/* contrib/recommend/recommend--1.0.sql */

-- complain if script is sourced in psql, rather than via ALTER EXTENSION
\echo Use "CREATE EXTENSION recommend" to load this file. \quit

CREATE OR REPLACE FUNCTION arraysml(anyarray, anyarray)
RETURNS	float4
AS	'MODULE_PATHNAME', 'arraysml'
LANGUAGE C STRICT IMMUTABLE;

CREATE OR REPLACE FUNCTION similarity(text, text)
RETURNS	float4
AS	'MODULE_PATHNAME', 'similarity'
LANGUAGE C STRICT IMMUTABLE;


CREATE FUNCTION twoint(integer,integer)
RETURNS integer
AS 'MODULE_PATHNAME', 'twoint'
LANGUAGE C STRICT IMMUTABLE;

CREATE FUNCTION two_float(float,float)
RETURNS float8
AS 'MODULE_PATHNAME', 'two_float'
LANGUAGE C STRICT IMMUTABLE;
