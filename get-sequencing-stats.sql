
SELECT 
	MAX(CASE process_udf_view.udfname WHEN 'Run ID' THEN process_udf_view.udfvalue ELSE '' END) AS "Run ID",
	project.name,
	artifact.name,
	MAX(CASE process_udf_view.udfname WHEN 'Read 1 Cycles' THEN process_udf_view.udfvalue ELSE '' END) AS "Read 1 Cycles",
	MAX(CASE process_udf_view.udfname WHEN 'Read 2 Cycles' THEN process_udf_view.udfvalue ELSE '' END) AS "Read 2 Cycles",
	MAX(CASE process_udf_view.udfname WHEN 'Index 1 Read Cycles' THEN process_udf_view.udfvalue ELSE '' END) AS "Index 1 Read Cycles",
	MAX(CASE process_udf_view.udfname WHEN 'Index 2 Read Cycles' THEN process_udf_view.udfvalue ELSE '' END) AS "Index 2 Read Cycles",
	MAX(CASE artifact_udf_view.udfname WHEN 'Yield PF (Gb) R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Yield PF (Gb) R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Yield PF (Gb) R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Yield PF (Gb) R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Avg Q Score R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Avg Q Score R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Avg Q Score R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Avg Q Score R2",
	MAX(CASE artifact_udf_view.udfname WHEN '% Bases >=Q30 R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Bases >=Q30 R1",
	MAX(CASE artifact_udf_view.udfname WHEN '% Bases >=Q30 R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Bases >=Q30 R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Cluster Density (K/mm^2) R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Cluster Density (K/mm^2) R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Cluster Density (K/mm^2) R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Cluster Density (K/mm^2) R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Clusters Raw R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Clusters Raw R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Clusters Raw R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Clusters Raw R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Clusters PF R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Clusters PF R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Clusters PF R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Clusters PF R2",
	MAX(CASE artifact_udf_view.udfname WHEN '%PF R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "%PF R1",
	MAX(CASE artifact_udf_view.udfname WHEN '%PF R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "%PF R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Intensity Cycle 1 R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Intensity Cycle 1 R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Intensity Cycle 1 R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Intensity Cycle 1 R2",
	MAX(CASE artifact_udf_view.udfname WHEN 'Intensity Cycle 20 R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Intensity Cycle 20 R1",
	MAX(CASE artifact_udf_view.udfname WHEN 'Intensity Cycle 20 R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "Intensity Cycle 20 R2",
	MAX(CASE artifact_udf_view.udfname WHEN '% Phasing R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Phasing R1",
	MAX(CASE artifact_udf_view.udfname WHEN '% Phasing R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Phasing R2",
	MAX(CASE artifact_udf_view.udfname WHEN '% Prephasing R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Prephasing R1",
	MAX(CASE artifact_udf_view.udfname WHEN '% Prephasing R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Prephasing R2",
	MAX(CASE artifact_udf_view.udfname WHEN '% Aligned R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Aligned R1",
	MAX(CASE artifact_udf_view.udfname WHEN '% Aligned R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Aligned R2",
	MAX(CASE artifact_udf_view.udfname WHEN '% Error Rate R1' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Error Rate R1",
	MAX(CASE artifact_udf_view.udfname WHEN '% Error Rate R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Error Rate R2",
	artifactstate.qcflag

FROM processtype, process, processiotracker, artifact, artifact_sample_map, sample, project, process_udf_view, artifact_udf_view, artifactstate
WHERE
	processtype.displayname='Illumina Sequencing (Illumina SBS) 5.0' AND
	--processtype.displayname='MiSeq Run (MiSeq) 5.0' AND--
	--processtype.displayname='NextSeq Run (NextSeq) 1.0' AND--
	process.typeid=processtype.typeid AND
	processiotracker.processid=process.processid AND
	artifact.artifactid=processiotracker.inputartifactid AND
	artifact_udf_view.artifactid=artifact.artifactid AND
	process_udf_view.processid=process.processid AND
	artifact_sample_map.artifactid=artifact.artifactid AND
	artifact_sample_map.processid=sample.processid AND
	project.projectid=sample.projectid AND
	artifactstate.artifactid=artifact.artifactid AND
	artifactstate.qcflag != 0

	GROUP BY process.processid, artifact.artifactid, artifact.name, project.name, artifact_sample_map.artifactid, artifactstate.qcflag
	ORDER BY "Run ID"


