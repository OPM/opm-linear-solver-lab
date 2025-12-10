.mode col
ATTACH 'bsr_test.db' as bsr;
--SELECT NAME FROM pragma_table_info("KERN");


SELECT NAME FROM BSR.SQLITE_MASTER;
CREATE TABLE STATS AS 
SELECT 
--Index
KernelName,
--gpu-id
--queue-id
--queue-index
--pid
--tid
count() as calls,
grd,
wgr,
lds,
scr,
arch_vgpr as vgpr,
accum_vgpr as agpr,
sgpr,
L2CacheHit as l2hits,
--LDSBankConflict, 
--wave_size, 64
--sig
--obj
--SUM(FETCH_SIZE)*1024 as FETCH_SIZE, -- *1024 = bytes
--SUM(WRITE_SIZE)*1024 as WRITE_SIZE,
--TCP_PENDING_STALL_CYCLES_sum as cycles_stalled_on_L2,
SUM(FETCH_SIZE)/1024 as READ_MB, -- *1024 = bytes
SUM(WRITE_SIZE)/1024 as WRITE_MB,
--DispatchNs
--BeginNs
--EndNs
--CompleteNs
--VALUUtilization as valutilization,
SUM(DurationNs) as DurationNs 
--FROM KERN 
FROM KERN WHERE KERNELNAME GLOB "*bsr*"
GROUP BY KERNELNAME;

SELECT 
trim(KernelName,'"') as kernel_name,
calls, 
PRINTF("%6.2f",READ_MB) as rd_mb,
PRINTF("%6.2f",WRITE_MB) as wt_mb,
DurationNs as DurationNs,
DurationNs/calls as duration_per_callNs,
--LDSBankConflict,
--valutilization,
l2hits,
--PRINTF("%6.2f",LDSBankConflict) as LDSBankConflict,
PRINTF("%6.2f",READ_MB/1024/(DurationNs*1e-9)) as read_bw_GBPS,
PRINTF("%6.2f",WRITE_MB/1024/(DurationNs*1e-9)) as write_bw_GBPS
FROM STATS;
--SELECT NAME FROM pragma_table_info("KERN");
DETACH bsr;

  
