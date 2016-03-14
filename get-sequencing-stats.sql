
SELECT 
	process_udf_view.udfvalue AS "Run ID",
	project.name,
	artifact.name,
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
	MAX(CASE artifact_udf_view.udfname WHEN '% Error Rate R2' THEN artifact_udf_view.udfvalue ELSE '' END) AS "% Error Rate R2"
FROM processtype, process, processiotracker, artifact, artifact_sample_map, sample, project, process_udf_view, artifact_udf_view
WHERE
	processtype.displayname='Illumina Sequencing (Illumina SBS) 5.0' AND
	process.typeid=processtype.typeid AND
	processiotracker.processid=process.processid AND
	artifact.artifactid=processiotracker.inputartifactid AND
	artifact_udf_view.artifactid=artifact.artifactid AND
	process_udf_view.processid=process.processid AND
	artifact_sample_map.artifactid=artifact.artifactid AND
	artifact_sample_map.processid=sample.processid AND
	project.projectid=sample.projectid AND
	process_udf_view.udfname='Run ID' AND

	0 < (SELECT COUNT(*) FROM project project2, sample sample2, artifact_sample_map asm2, artifact artifact2, processiotracker processiotracker2
		WHERE
		processiotracker2.processid=process.processid AND
		processiotracker2.inputartifactid=artifact2.artifactid AND
		asm2.artifactid=artifact2.artifactid AND
		asm2.processid=sample2.processid AND
		project2.projectid=sample2.projectid AND
		UPPER(project2.name) LIKE 'DIAG%')

	GROUP BY process.processid, artifact.artifactid, process_udf_view.udfvalue, artifact.name, project.name
	ORDER BY process_udf_view.udfvalue

