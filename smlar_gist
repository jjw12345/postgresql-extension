#include "smlar.h"

#include "fmgr.h"
#include "access/gist.h"
#include "access/skey.h"
#include "access/tuptoaster.h"
#include "utils/memutils.h"

typedef struct SmlSign {
	int32	vl_len_; /* varlena header (do not touch directly!) */
	int32	flag:8,
			size:24;
	int32	maxrepeat;
	char	data[1];
} SmlSign;

#define SMLSIGNHDRSZ	(offsetof(SmlSign, data))

#define BITBYTE 8
#define SIGLENINT  61
#define SIGLEN  ( sizeof(int)*SIGLENINT )
#define SIGLENBIT (SIGLEN*BITBYTE - 1)  /* see makesign */
typedef char BITVEC[SIGLEN];
typedef char *BITVECP;
#define LOOPBYTE \
		for(i=0;i<SIGLEN;i++)

#define GETBYTE(x,i) ( *( (BITVECP)(x) + (int)( (i) / BITBYTE ) ) )
#define GETBITBYTE(x,i) ( ((char)(x)) >> i & 0x01 )
#define CLRBIT(x,i)   GETBYTE(x,i) &= ~( 0x01 << ( (i) % BITBYTE ) )
#define SETBIT(x,i)   GETBYTE(x,i) |=  ( 0x01 << ( (i) % BITBYTE ) )
#define GETBIT(x,i) ( (GETBYTE(x,i) >> ( (i) % BITBYTE )) & 0x01 )

#define HASHVAL(val) (((unsigned int)(val)) % SIGLENBIT)
#define HASH(sign, val) SETBIT((sign), HASHVAL(val))

#define ARRKEY          0x01
#define SIGNKEY         0x02
#define ALLISTRUE       0x04

#define ISARRKEY(x) ( ((SmlSign*)x)->flag & ARRKEY )
#define ISSIGNKEY(x)    ( ((SmlSign*)x)->flag & SIGNKEY )
#define ISALLTRUE(x)    ( ((SmlSign*)x)->flag & ALLISTRUE )

#define CALCGTSIZE(flag, len) ( SMLSIGNHDRSZ + ( ( (flag) & ARRKEY ) ? ((len)*sizeof(uint32)) : (((flag) & ALLISTRUE) ? 0 : SIGLEN) ) )
#define GETSIGN(x)      ( (BITVECP)( (char*)x+SMLSIGNHDRSZ ) )
#define GETARR(x)       ( (uint32*)( (char*)x+SMLSIGNHDRSZ ) )

#define GETENTRY(vec,pos) ((SmlSign *) DatumGetPointer((vec)->vector[(pos)].key))

/*
 * Fake IO
 */

/*
 * Compress/decompress
 */

/* Number of one-bits in an unsigned byte */
static const uint8 number_of_ones[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

static int
compareint(const void *va, const void *vb)
{
    uint32        a = *((uint32 *) va);
	uint32        b = *((uint32 *) vb);
				 
	if (a == b)
		return 0;
	return (a > b) ? 1 : -1;
}

/*
 * Removes duplicates from an array of int32. 'l' is
 * size of the input array. Returns the new size of the array.
 */
static int
uniqueint(uint32 *a, int32 l, int32 *max)
{
	uint32      *ptr,
			   *res;
	int32	   cnt = 0;

	*max = 1;

	if (l <= 1)
		return l;

	ptr = res = a;

	qsort((void *) a, l, sizeof(uint32), compareint);

	while (ptr - a < l)
		if (*ptr != *res)
		{
			cnt = 1;
			*(++res) = *ptr++;
		}
		else
		{
			cnt++;
			if ( cnt > *max )
				*max = cnt;
			ptr++;
		}

	if ( cnt > *max )
		*max = cnt;

	return res + 1 - a;
}

SmlSign*
Array2HashedArray(ProcTypeInfo info, ArrayType *a)
{
	SimpleArray *s = Array2SimpleArray(info, a);
	SmlSign		*sign;
	int32		len, i;
	uint32		*ptr;

	len = CALCGTSIZE( ARRKEY, s->nelems );

	getFmgrInfoHash(s->info);
	if (s->info->tupDesc)
		elog(ERROR, "GiST  doesn't support composite (weighted) type");

	sign = palloc( len );
	sign->flag = ARRKEY;
	sign->size = s->nelems;

	ptr = GETARR(sign);
	for(i=0;i<s->nelems;i++)
		ptr[i] = DatumGetUInt32( FunctionCall1( &s->info->hashFunc, s->elems[i] ) );
				
	/*
	 * there is a collision of hash-function; len is always equal or less than
	 * s->nelems
	 */
	sign->size = uniqueint( GETARR(sign), sign->size, &sign->maxrepeat );
	len = CALCGTSIZE( ARRKEY, sign->size );
	SET_VARSIZE(sign, len);

	return sign;
}

#define WISH_F(a,b,c) (double)( -(double)(((a)-(b))*((a)-(b))*((a)-(b)))*(c) )

typedef struct
{
	OffsetNumber pos;
	int32        cost;
} SPLITCOST;
