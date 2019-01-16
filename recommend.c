/*-------------------------------------------------------------------------
 *
 * tcn.c
 *	  triggered change notification support for PostgreSQL
 *
 * Portions Copyright (c) 2011-2018, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 *
 * IDENTIFICATION
 *	  contrib/recommend/recommend.c
 *
 *-------------------------------------------------------------------------
 */

#include "smlar.h"

#include "fmgr.h"
#include "access/genam.h"
#include "access/heapam.h"
#include "access/htup_details.h"
#include "access/nbtree.h"
#include "catalog/indexing.h"
#include "catalog/pg_am.h"
#include "catalog/pg_amproc.h"
#include "catalog/pg_cast.h"
#include "catalog/pg_opclass.h"
#include "catalog/pg_type.h"
#include "executor/spi.h"
#include "utils/catcache.h"
#include "utils/fmgroids.h"
#include "utils/lsyscache.h"
#include "utils/memutils.h"
#include "utils/tqual.h"
#include "utils/syscache.h"
#include "utils/typcache.h"
#include "trgm.h"
#include<stdio.h>
#include<malloc.h>
#include<math.h>


#if (PG_VERSION_NUM >= 90400)
#define SNAPSHOT NULL
#else
#define SNAPSHOT SnapshotNow
#endif

PG_MODULE_MAGIC;


static Oid
getDefaultOpclass(Oid amoid, Oid typid)
{
	ScanKeyData	skey;
	SysScanDesc	scan;
	HeapTuple	tuple;
	Relation	heapRel;
	Oid			opclassOid = InvalidOid;

	heapRel = heap_open(OperatorClassRelationId, AccessShareLock);

	ScanKeyInit(&skey,
				Anum_pg_opclass_opcmethod,
				BTEqualStrategyNumber,	F_OIDEQ,
				ObjectIdGetDatum(amoid));

	scan = systable_beginscan(heapRel,
								OpclassAmNameNspIndexId, true,
								SNAPSHOT, 1, &skey);

	while (HeapTupleIsValid((tuple = systable_getnext(scan))))
	{
		Form_pg_opclass opclass = (Form_pg_opclass)GETSTRUCT(tuple);

		if ( opclass->opcintype == typid && opclass->opcdefault )
		{
			if ( OidIsValid(opclassOid) )
				elog(ERROR, "Ambiguous opclass for type %u (access method %u)", typid, amoid); 
			opclassOid = HeapTupleGetOid(tuple);
		}
	}

	systable_endscan(scan);
	heap_close(heapRel, AccessShareLock);

	return opclassOid;
}

static Oid
getAMProc(Oid amoid, Oid typid)
{
	Oid		opclassOid = getDefaultOpclass(amoid, typid);
	Oid		procOid = InvalidOid;
	Oid		opfamilyOid;
	ScanKeyData	skey[4];
	SysScanDesc	scan;
	HeapTuple	tuple;
	Relation	heapRel;

	if ( !OidIsValid(opclassOid) )
	{
		typid = getBaseType(typid);
		opclassOid = getDefaultOpclass(amoid, typid);
	}

	if ( !OidIsValid(opclassOid) )
	{
		CatCList	*catlist;
		int			i;

		/*
		 * Search binary-coercible type
		 */
		catlist = SearchSysCacheList(CASTSOURCETARGET, 1,
										ObjectIdGetDatum(typid),
										0, 0, 0);

		for (i = 0; i < catlist->n_members; i++)
		{
			HeapTuple		tuple = &catlist->members[i]->tuple;
			Form_pg_cast	castForm = (Form_pg_cast)GETSTRUCT(tuple);

			if ( castForm->castfunc == InvalidOid && castForm->castcontext == COERCION_CODE_IMPLICIT )
			{
				typid = castForm->casttarget;
				opclassOid = getDefaultOpclass(amoid, typid);
				if( OidIsValid(opclassOid) )
					break;
			}
		}

		ReleaseSysCacheList(catlist);
	}

	if ( !OidIsValid(opclassOid) )
		return InvalidOid;

	opfamilyOid = get_opclass_family(opclassOid);

	heapRel = heap_open(AccessMethodProcedureRelationId, AccessShareLock);
	ScanKeyInit(&skey[0],
				Anum_pg_amproc_amprocfamily,
				BTEqualStrategyNumber, F_OIDEQ,
				ObjectIdGetDatum(opfamilyOid));
	ScanKeyInit(&skey[1],
				Anum_pg_amproc_amproclefttype,
				BTEqualStrategyNumber, F_OIDEQ,
				ObjectIdGetDatum(typid));
	ScanKeyInit(&skey[2],
				Anum_pg_amproc_amprocrighttype,
				BTEqualStrategyNumber, F_OIDEQ,
				ObjectIdGetDatum(typid));
#if PG_VERSION_NUM >= 90200
	ScanKeyInit(&skey[3],
				Anum_pg_amproc_amprocnum,
				BTEqualStrategyNumber, F_OIDEQ,
				Int32GetDatum(BTORDER_PROC));
#endif

	scan = systable_beginscan(heapRel, AccessMethodProcedureIndexId, true,
								SNAPSHOT,
#if PG_VERSION_NUM >= 90200
								4,
#else
								3,
#endif
								skey);
	while (HeapTupleIsValid(tuple = systable_getnext(scan)))
	{
		Form_pg_amproc amprocform = (Form_pg_amproc) GETSTRUCT(tuple);

		switch(amoid)
		{
			case BTREE_AM_OID:
			case HASH_AM_OID:
				if ( OidIsValid(procOid) )
					elog(ERROR,"Ambiguous support function for type %u (opclass %u)", typid, opfamilyOid);
				procOid = amprocform->amproc;
				break;
			default:
				elog(ERROR,"Unsupported access method");
		}
	}

	systable_endscan(scan);
	heap_close(heapRel, AccessShareLock);

	return procOid;
}

static ProcTypeInfo *cacheProcs = NULL;
static int nCacheProcs = 0;

static ProcTypeInfo
fillProcs(Oid typid)
{
	ProcTypeInfo	info = malloc(sizeof(ProcTypeInfoData));

	if (!info)
		elog(ERROR, "Can't allocate %u memory", (uint32)sizeof(ProcTypeInfoData));

	info->typid = typid;
	info->typtype = get_typtype(typid);

	if (info->typtype == 'c')
	{
		/* composite type */
		TupleDesc		tupdesc;
		MemoryContext	oldcontext;

		tupdesc = lookup_rowtype_tupdesc(typid, -1);

		if (tupdesc->natts != 2)
			elog(ERROR,"Composite type has wrong number of fields");
		if (tupdesc->attrs[1]->atttypid != FLOAT4OID)
			elog(ERROR,"Second field of composite type is not float4");

		oldcontext = MemoryContextSwitchTo(TopMemoryContext);
		info->tupDesc = CreateTupleDescCopyConstr(tupdesc);
		MemoryContextSwitchTo(oldcontext);

		ReleaseTupleDesc(tupdesc);

		info->cmpFuncOid = getAMProc(BTREE_AM_OID, info->tupDesc->attrs[0]->atttypid);
		info->hashFuncOid = getAMProc(HASH_AM_OID, info->tupDesc->attrs[0]->atttypid);
	}
	else
	{
		info->tupDesc = NULL;

		/* plain type */
		info->cmpFuncOid = getAMProc(BTREE_AM_OID, typid);
		info->hashFuncOid = getAMProc(HASH_AM_OID, typid);
	}

	get_typlenbyvalalign(typid, &info->typlen, &info->typbyval, &info->typalign);
	info->hashFuncInited = info->cmpFuncInited = false;


	return info;
}

void
getFmgrInfoCmp(ProcTypeInfo info)
{
	if ( info->cmpFuncInited == false )
	{
		if ( !OidIsValid(info->cmpFuncOid) )
			elog(ERROR, "Could not find cmp function for type %u", info->typid);

		fmgr_info_cxt( info->cmpFuncOid, &info->cmpFunc, TopMemoryContext );
		info->cmpFuncInited = true;
	}
}

void
getFmgrInfoHash(ProcTypeInfo info)
{
	if ( info->hashFuncInited == false )
	{
		if ( !OidIsValid(info->hashFuncOid) )
			elog(ERROR, "Could not find hash function for type %u", info->typid);

		fmgr_info_cxt( info->hashFuncOid, &info->hashFunc, TopMemoryContext );
		info->hashFuncInited = true;
	}
}

static int
cmpProcTypeInfo(const void *a, const void *b)
{
	ProcTypeInfo av = *(ProcTypeInfo*)a;
	ProcTypeInfo bv = *(ProcTypeInfo*)b;

	Assert( av->typid != bv->typid );

	return ( av->typid > bv->typid ) ? 1 : -1;
}

ProcTypeInfo
findProcs(Oid typid)
{
	ProcTypeInfo	info = NULL;

	if ( nCacheProcs == 1 )
	{
		if ( cacheProcs[0]->typid == typid )
		{
			/*cacheProcs[0]->hashFuncInited = cacheProcs[0]->cmpFuncInited = false;*/
			return cacheProcs[0];
		}
	}
	else if ( nCacheProcs > 1 )
	{
		ProcTypeInfo	*StopMiddle;
		ProcTypeInfo	*StopLow = cacheProcs,
						*StopHigh = cacheProcs + nCacheProcs;

		while (StopLow < StopHigh) {
			StopMiddle = StopLow + ((StopHigh - StopLow) >> 1);
			info = *StopMiddle;

			if ( info->typid == typid )
			{
				/* info->hashFuncInited = info->cmpFuncInited = false; */
				return info;
			}
			else if ( info->typid < typid )
				StopLow = StopMiddle + 1;
			else
				StopHigh = StopMiddle;
		}

		/* not found */
	} 

	info = fillProcs(typid);
	if ( nCacheProcs == 0 )
	{
		cacheProcs = malloc(sizeof(ProcTypeInfo));

		if (!cacheProcs)
			elog(ERROR, "Can't allocate %u memory", (uint32)sizeof(ProcTypeInfo));
		else
		{
			nCacheProcs = 1;
			cacheProcs[0] = info;
		}
	}
	else
	{
		ProcTypeInfo	*cacheProcsTmp = realloc(cacheProcs, (nCacheProcs+1) * sizeof(ProcTypeInfo));

		if (!cacheProcsTmp)
			elog(ERROR, "Can't allocate %u memory", (uint32)sizeof(ProcTypeInfo) * (nCacheProcs+1));
		else
		{
			cacheProcs = cacheProcsTmp;
			cacheProcs[ nCacheProcs ] = info;
			nCacheProcs++;
			qsort(cacheProcs, nCacheProcs, sizeof(ProcTypeInfo), cmpProcTypeInfo);
		}
	}

	/* info->hashFuncInited = info->cmpFuncInited = false; */

	return info;
}


/*
 * WARNING. Array2SimpleArray* doesn't copy Datum!
 */
SimpleArray * 
Array2SimpleArray(ProcTypeInfo info, ArrayType *a)
{
	SimpleArray	*s = palloc(sizeof(SimpleArray));

	CHECKARRVALID(a);

	if ( info == NULL )
		info = findProcs(ARR_ELEMTYPE(a));

	s->info = info;
	s->df = NULL;
	s->hash = NULL;

	deconstruct_array(a, info->typid,
						info->typlen, info->typbyval, info->typalign,
						&s->elems, NULL, &s->nelems);

	return s;
}

//构建类型 
static Datum
deconstructCompositeType(ProcTypeInfo info, Datum in, double *weight)
{
	HeapTupleHeader	rec = DatumGetHeapTupleHeader(in);
	HeapTupleData	tuple;
	Datum			values[2];
	bool			nulls[2];

	/* Build a temporary HeapTuple control structure */
	tuple.t_len = HeapTupleHeaderGetDatumLength(rec);   
	ItemPointerSetInvalid(&(tuple.t_self));    
	tuple.t_tableOid = InvalidOid;
	tuple.t_data = rec;

	heap_deform_tuple(&tuple, info->tupDesc, values, nulls);    
	if (nulls[0] || nulls[1])
		elog(ERROR, "Both fields in composite type could not be NULL");

	if (weight)
		*weight = DatumGetFloat4(values[1]);   
	return values[0];
}





static int
cmpArrayElem(const void *a, const void *b, void *arg)
{
	ProcTypeInfo	info = (ProcTypeInfo)arg;  

	if (info->tupDesc)
		/* composite type */
		return DatumGetInt32( FCall2( &info->cmpFunc,
						deconstructCompositeType(info, *(Datum*)a, NULL),
						deconstructCompositeType(info, *(Datum*)b, NULL) ) );   
	return DatumGetInt32( FCall2( &info->cmpFunc,
							*(Datum*)a, *(Datum*)b ) );    
}


SimpleArray *
Array2SimpleArrayS(ProcTypeInfo info, ArrayType *a)
{
	SimpleArray	*s = Array2SimpleArray(info, a);

	if ( s->nelems > 1 )
	{
		getFmgrInfoCmp(s->info);

		qsort_arg(s->elems, s->nelems, sizeof(Datum), cmpArrayElem, s->info);
	}

	return s;
}

typedef struct cmpArrayElemData {
	ProcTypeInfo	info;
	bool			hasDuplicate;

} cmpArrayElemData;

static int
cmpArrayElemArg(const void *a, const void *b, void *arg)
{
    cmpArrayElemData    *data = (cmpArrayElemData*)arg;
	int res;

	if (data->info->tupDesc)
		res =  DatumGetInt32( FCall2( &data->info->cmpFunc,
											deconstructCompositeType(data->info, *(Datum*)a, NULL),
											deconstructCompositeType(data->info, *(Datum*)b, NULL) ) );
	else
		res = DatumGetInt32( FCall2( &data->info->cmpFunc,
								*(Datum*)a, *(Datum*)b ) );

	if ( res == 0 )
		data->hasDuplicate = true;

	return res;
}

//新加
SimpleArray *
Array2SimpleArrayU(ProcTypeInfo info, ArrayType *a, void *cache)
{
	SimpleArray	*s = Array2SimpleArray(info, a);
	StatElem	*stat = NULL;

	if ( s->nelems > 0 && cache )
	{
		s->df = palloc(sizeof(double) * s->nelems);
		s->df[0] = 1.0; /* init */
	}

	if ( s->nelems > 1 )
	{
		cmpArrayElemData	data;
		int 				i;

		getFmgrInfoCmp(s->info);
		data.info = s->info;
		data.hasDuplicate = false;

		qsort_arg(s->elems, s->nelems, sizeof(Datum), cmpArrayElemArg, &data);

		if ( data.hasDuplicate )
		{
			Datum	*tmp,
					*dr,
					*data;
			int		num = s->nelems,
					cmp;

			data = tmp = dr = s->elems;

			while (tmp - data < num)
			{
				cmp = (tmp == dr) ? 0 : cmpArrayElem(tmp, dr, s->info);
				if ( cmp != 0 )
				{
					*(++dr) = *tmp++;
					if ( cache ) 
						s->df[ dr - data ] = 1.0;
				}
				else
				{
					if ( cache )
						s->df[ dr - data ] += 1.0;
					tmp++;
				}
			}

			s->nelems = dr + 1 - s->elems;

			if ( cache )
			{
				int tfm = getTFMethod();

				for(i=0;i<s->nelems;i++)
				{
					stat = fingArrayStat(cache, s->info->typid, s->elems[i], stat);
					if ( stat )
					{
						switch(tfm)
						{
							case TF_LOG:
								s->df[i] = (1.0 + log( s->df[i] ));
							case TF_N:
								s->df[i] *= stat->idf;
								break;
							case TF_CONST:
								s->df[i] = stat->idf;
								break;
							default:
								elog(ERROR,"Unknown TF method: %d", tfm);
						}	
					}
					else
					{
						s->df[i] = 0.0; /* unknown word */
					}
				}	
			}
		}
		else if ( cache )
		{
			for(i=0;i<s->nelems;i++)
			{
				stat = fingArrayStat(cache, s->info->typid, s->elems[i], stat);
				if ( stat )
					s->df[i] = stat->idf;
				else
					s->df[i] = 0.0;
			}
		}
	}
	else if (s->nelems > 0 && cache) 
	{
		stat = fingArrayStat(cache, s->info->typid, s->elems[0], stat);
		if ( stat )
			s->df[0] = stat->idf;
		else
			s->df[0] = 0.0;
	}

	return s;
}
 


//数组交叉个数的统计 
static int
numOfIntersect(SimpleArray *a, SimpleArray *b)    //数组交集个数统计 
{
	int				cnt = 0,
					cmp;
	Datum			*aptr = a->elems,    
					*bptr = b->elems;    
	ProcTypeInfo	info = a->info;      

	Assert( a->info->typid == b->info->typid );    

	getFmgrInfoCmp(info);   

	while( aptr - a->elems < a->nelems && bptr - b->elems < b->nelems )
	{
		cmp = cmpArrayElem(aptr, bptr, info);   
		if ( cmp < 0 )
			aptr++;
		else if ( cmp > 0 )
			bptr++;
		else
		{
			cnt++;  //记数 
			aptr++;
			bptr++;
		}
	}

	return cnt;
}
//tfidsml函数 
static double
TFIDFSml(SimpleArray *a, SimpleArray *b)
{
	int				cmp;
	Datum			*aptr = a->elems,
					*bptr = b->elems;
	ProcTypeInfo	info = a->info;
	double			res = 0.0;
	double			suma = 0.0, sumb = 0.0;

	Assert( a->info->typid == b->info->typid );
	Assert( a->df );
	Assert( b->df );

	getFmgrInfoCmp(info);

	while( aptr - a->elems < a->nelems && bptr - b->elems < b->nelems )
	{
		cmp = cmpArrayElem(aptr, bptr, info);
		if ( cmp < 0 )
		{
			suma += a->df[ aptr - a->elems ] * a->df[ aptr - a->elems ];
			aptr++;
		}
		else if ( cmp > 0 )
		{
			sumb += b->df[ bptr - b->elems ] * b->df[ bptr - b->elems ];
			bptr++;
		}
		else
		{
			res += a->df[ aptr - a->elems ] * b->df[ bptr - b->elems ];
			suma += a->df[ aptr - a->elems ] * a->df[ aptr - a->elems ];
			sumb += b->df[ bptr - b->elems ] * b->df[ bptr - b->elems ];
			aptr++;
			bptr++;
		}
	}

	/*
	 * Compute last elements
	 */
	while( aptr - a->elems < a->nelems )
	{
		suma += a->df[ aptr - a->elems ] * a->df[ aptr - a->elems ];
		aptr++;
	}

	while( bptr - b->elems < b->nelems )
	{
		sumb += b->df[ bptr - b->elems ] * b->df[ bptr - b->elems ];
		bptr++;
	}

	if ( suma > 0.0 && sumb > 0.0 )
		res = res / sqrt( suma * sumb );
	else
		res = 0.0;

	return res;
}


//数组相似度  计算(?)
PG_FUNCTION_INFO_V1(arraysml);
Datum	arraysml(PG_FUNCTION_ARGS);
Datum
arraysml(PG_FUNCTION_ARGS)
{
	ArrayType		*a, *b;    //数组类型 
	SimpleArray		*sa, *sb;   //简单数组 

	fcinfo->flinfo->fn_extra = SearchArrayCache(
							fcinfo->flinfo->fn_extra,
							fcinfo->flinfo->fn_mcxt,
							PG_GETARG_DATUM(0), &a, &sa, NULL);     //属于smlar.h下的函数 
	fcinfo->flinfo->fn_extra = SearchArrayCache(
							fcinfo->flinfo->fn_extra,
							fcinfo->flinfo->fn_mcxt,
							PG_GETARG_DATUM(1), &b, &sb, NULL);

	if ( ARR_ELEMTYPE(a) != ARR_ELEMTYPE(b) )    //这个函数   应该属于头文件中的一个函数 
		elog(ERROR,"Arguments array are not the same type!");

	if (ARRISVOID(a) || ARRISVOID(b))
		 PG_RETURN_FLOAT4(0.0);

	switch(getSmlType())   
	{
		case ST_TFIDF:
			PG_RETURN_FLOAT4( TFIDFSml(sa, sb) );   
			break;
		case ST_COSINE:
			{
				int				cnt;
				double			power;

				power = ((double)(sa->nelems)) * ((double)(sb->nelems));   
				cnt = numOfIntersect(sa, sb);      

				PG_RETURN_FLOAT4(  ((double)cnt) / sqrt( power ) );    
			}
			break;
		case ST_OVERLAP:
			{
				float4 res = (float4)numOfIntersect(sa, sb);      

				PG_RETURN_FLOAT4(res);    
			}
			break;
		default:
			elog(ERROR,"Unsupported formula type of similarity");
	}

	PG_RETURN_FLOAT4(0.0); /* keep compiler quiet */
}


//字符串相似度  计算。 

typedef struct vect {
	int32		vl_len_;		/* varlena header (do not touch directly!) */
	uint8		flag;
	int32		data[28];
} vect;  

//向量计算 

float4 cossml(vect *trg1, vect *trg2)
{
	float4 g=0.0;
	float4 s1,s2;
	int i=0;
	s1=0;
	s2=0;
	for(i=0;i<28;i++)
	{
		s1=s1+ trg1->data[i]*trg1->data[i];
		s2=s2+trg2->data[i]*trg2->data[i];
		g=g+trg1->data[i]*trg2->data[i];
	}
	g=g/((float4)sqrt(s1)*(float4)sqrt(s2));
	return g;
}

//生成向量型。 
vect * generate_vect(char *a,int slen)   
	vect *c;
	int i;
	int m=0;
	
	c=(vect *)palloc(TRGMHDRSIZE+98);
	c->flag = ARRKEY; 
	
	for(i=0;i<28;i++)
	{
		c->data[i]=0;
	} 
	
	for(i=0;i<slen;i++)    
	{
		m=a[i]-'a';
		c->data[m]++;
	}
	c->vl_len_=slen;
	
	return c;
	
} 


PG_FUNCTION_INFO_V1(similarity);
Datum	similarity(PG_FUNCTION_ARGS);

Datum
similarity(PG_FUNCTION_ARGS)
{
	text	   *in1 = PG_GETARG_TEXT_PP(0);  
	text	   *in2 = PG_GETARG_TEXT_PP(1);  
	vect	   *trg1,
			   *trg2;                          
	float4		res;     

	//trg1 = generate_trgm(VARDATA_ANY(in1), VARSIZE_ANY_EXHDR(in1));
	//trg2 = generate_trgm(VARDATA_ANY(in2), VARSIZE_ANY_EXHDR(in2));
	trg1=generate_vect(VARDATA_ANY(in1), VARSIZE_ANY_EXHDR(in1));  //vardata_any(in1)
	trg2=generate_vect(VARDATA_ANY(in2), VARSIZE_ANY_EXHDR(in2));  

	res = cossml(trg1, trg2);    
	pfree(trg1);  //要加头文件. 
	pfree(trg2);
	//pfree(trg1);   //释放内存 
	//pfree(trg2);
	//PG_FREE_IF_COPY(in1, 0);
	//PG_FREE_IF_COPY(in2, 1);

	PG_RETURN_FLOAT4(res);       
}







//函数 int 
PG_FUNCTION_INFO_V1(twoint);  
Datum twoint(PG_FUNCTION_ARGS);   
Datum twoint(PG_FUNCTION_ARGS)
{
	int sum ,a,b;
	
	a=PG_GETARG_INT32(0);   
	b=PG_GETARG_INT32(1);
    if(a>b)
    {
    	sum = a - b;
	}
	else
	{
		sum=b-a;
	}
	
	PG_RETURN_INT32(sum);      
}

//函数float 
PG_FUNCTION_INFO_V1(two_float);
Datum	two_float(PG_FUNCTION_ARGS);
Datum
two_float(PG_FUNCTION_ARGS){
	float8 sum,a,b;
	a=PG_GETARG_FLOAT8(0);     
	b=PG_GETARG_FLOAT8(1);     
	if(a>b)
	{
		sum=a-b;
	}
	else
	{
		sum=b-a;
	}
	PG_RETURN_FLOAT8(sum);
}
