-- ============================================================
-- PharmacoDB — RPC Functions for Dashboard
-- Accessible via Supabase REST API (public schema)
-- ============================================================

-- Grant read access
GRANT USAGE ON SCHEMA pharmaco_db TO anon, authenticated;
GRANT SELECT ON ALL TABLES IN SCHEMA pharmaco_db TO anon, authenticated;
ALTER DEFAULT PRIVILEGES IN SCHEMA pharmaco_db GRANT SELECT ON TABLES TO anon, authenticated;

-- 1. Overview KPIs
CREATE OR REPLACE FUNCTION public.pharmaco_overview()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_build_object(
    'compounds', (SELECT COUNT(*) FROM pharmaco_db.compounds),
    'targets', (SELECT COUNT(*) FROM pharmaco_db.targets),
    'assays', (SELECT COUNT(*) FROM pharmaco_db.assays),
    'bioactivities', (SELECT COUNT(*) FROM pharmaco_db.bioactivities),
    'drug_mechanisms', (SELECT COUNT(*) FROM pharmaco_db.drug_mechanisms),
    'drug_indications', (SELECT COUNT(*) FROM pharmaco_db.drug_indications),
    'approved_drugs', (SELECT COUNT(*) FROM pharmaco_db.compounds WHERE max_phase >= 4),
    'druggable_targets', (SELECT COUNT(*) FROM pharmaco_db.targets WHERE is_druggable = TRUE),
    'targets_with_activity', (SELECT COUNT(DISTINCT target_id) FROM pharmaco_db.bioactivities WHERE target_id IS NOT NULL),
    'compounds_with_activity', (SELECT COUNT(DISTINCT compound_id) FROM pharmaco_db.bioactivities),
    'avg_pchembl', (SELECT ROUND(AVG(pchembl_value)::numeric, 2) FROM pharmaco_db.bioactivities WHERE pchembl_value IS NOT NULL),
    'db_size', (SELECT pg_size_pretty(pg_database_size('postgres')))
  ) INTO result;
  RETURN result;
END;
$$;

-- 2. Molecular weight distribution
CREATE OR REPLACE FUNCTION public.pharmaco_mw_distribution()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      width_bucket(molecular_weight, 100, 900, 16) AS bin,
      ROUND((100 + (width_bucket(molecular_weight, 100, 900, 16) - 1) * 50)::numeric) AS mw_low,
      ROUND((100 + width_bucket(molecular_weight, 100, 900, 16) * 50)::numeric) AS mw_high,
      COUNT(*) AS count
    FROM pharmaco_db.compounds
    WHERE molecular_weight IS NOT NULL
    GROUP BY bin ORDER BY bin
  ) t;
  RETURN result;
END;
$$;

-- 3. LogP distribution
CREATE OR REPLACE FUNCTION public.pharmaco_logp_distribution()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      width_bucket(alogp, -5, 10, 15) AS bin,
      ROUND((-5 + (width_bucket(alogp, -5, 10, 15) - 1))::numeric, 1) AS logp_low,
      ROUND((-5 + width_bucket(alogp, -5, 10, 15))::numeric, 1) AS logp_high,
      COUNT(*) AS count
    FROM pharmaco_db.compounds
    WHERE alogp IS NOT NULL
    GROUP BY bin ORDER BY bin
  ) t;
  RETURN result;
END;
$$;

-- 4. Target families breakdown
CREATE OR REPLACE FUNCTION public.pharmaco_target_families()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      COALESCE(protein_family, 'Unknown') AS family,
      COUNT(*) AS count
    FROM pharmaco_db.targets
    GROUP BY protein_family
    ORDER BY count DESC
    LIMIT 20
  ) t;
  RETURN result;
END;
$$;

-- 5. Activity type breakdown
CREATE OR REPLACE FUNCTION public.pharmaco_activity_types()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      activity_type AS type,
      COUNT(*) AS count,
      ROUND(AVG(pchembl_value)::numeric, 2) AS avg_pchembl
    FROM pharmaco_db.bioactivities
    WHERE pchembl_value IS NOT NULL
    GROUP BY activity_type
    ORDER BY count DESC
    LIMIT 15
  ) t;
  RETURN result;
END;
$$;

-- 6. pChEMBL distribution
CREATE OR REPLACE FUNCTION public.pharmaco_pchembl_distribution()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      width_bucket(pchembl_value, 3, 12, 18) AS bin,
      ROUND((3 + (width_bucket(pchembl_value, 3, 12, 18) - 1) * 0.5)::numeric, 1) AS pchembl_low,
      ROUND((3 + width_bucket(pchembl_value, 3, 12, 18) * 0.5)::numeric, 1) AS pchembl_high,
      COUNT(*) AS count
    FROM pharmaco_db.bioactivities
    WHERE pchembl_value IS NOT NULL
    GROUP BY bin ORDER BY bin
  ) t;
  RETURN result;
END;
$$;

-- 7. Clinical phase funnel
CREATE OR REPLACE FUNCTION public.pharmaco_clinical_phases()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      max_phase AS phase,
      CASE max_phase
        WHEN 0 THEN 'Preclinical'
        WHEN 1 THEN 'Phase I'
        WHEN 2 THEN 'Phase II'
        WHEN 3 THEN 'Phase III'
        WHEN 4 THEN 'Approved'
        ELSE 'Unknown'
      END AS label,
      COUNT(*) AS count
    FROM pharmaco_db.compounds
    WHERE max_phase IS NOT NULL
    GROUP BY max_phase
    ORDER BY max_phase
  ) t;
  RETURN result;
END;
$$;

-- 8. Chemical space (MW vs logP, sampled)
CREATE OR REPLACE FUNCTION public.pharmaco_chemical_space(p_limit int DEFAULT 2000)
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      molecular_weight AS mw,
      alogp AS logp,
      COALESCE(qed_weighted, 0) AS qed,
      COALESCE(max_phase, 0) AS phase,
      num_ro5_violations AS ro5v
    FROM pharmaco_db.compounds
    WHERE molecular_weight IS NOT NULL AND alogp IS NOT NULL
    ORDER BY RANDOM()
    LIMIT p_limit
  ) t;
  RETURN result;
END;
$$;

-- 9. Drug-likeness summary
CREATE OR REPLACE FUNCTION public.pharmaco_druglikeness()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_build_object(
    'ro5_distribution', (
      SELECT json_agg(row_to_json(t)) FROM (
        SELECT num_ro5_violations AS violations, COUNT(*) AS count
        FROM pharmaco_db.compounds
        WHERE num_ro5_violations IS NOT NULL
        GROUP BY num_ro5_violations ORDER BY num_ro5_violations
      ) t
    ),
    'qed_distribution', (
      SELECT json_agg(row_to_json(t)) FROM (
        SELECT
          width_bucket(qed_weighted, 0, 1, 10) AS bin,
          ROUND((width_bucket(qed_weighted, 0, 1, 10) - 1)::numeric * 0.1, 1) AS qed_low,
          ROUND(width_bucket(qed_weighted, 0, 1, 10)::numeric * 0.1, 1) AS qed_high,
          COUNT(*) AS count
        FROM pharmaco_db.compounds
        WHERE qed_weighted IS NOT NULL
        GROUP BY bin ORDER BY bin
      ) t
    ),
    'avg_qed', (SELECT ROUND(AVG(qed_weighted)::numeric, 3) FROM pharmaco_db.compounds WHERE qed_weighted IS NOT NULL),
    'pct_drug_like', (SELECT ROUND(100.0 * COUNT(*) FILTER (WHERE num_ro5_violations <= 1) / NULLIF(COUNT(*), 0), 1) FROM pharmaco_db.compounds WHERE num_ro5_violations IS NOT NULL)
  ) INTO result;
  RETURN result;
END;
$$;

-- 10. Mechanism of action types
CREATE OR REPLACE FUNCTION public.pharmaco_mechanisms()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      COALESCE(action_type, 'Unknown') AS action,
      COUNT(*) AS count
    FROM pharmaco_db.drug_mechanisms
    GROUP BY action_type
    ORDER BY count DESC
    LIMIT 15
  ) t;
  RETURN result;
END;
$$;

-- 11. Target types
CREATE OR REPLACE FUNCTION public.pharmaco_target_types()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT
      COALESCE(target_type, 'Unknown') AS type,
      COUNT(*) AS count
    FROM pharmaco_db.targets
    GROUP BY target_type
    ORDER BY count DESC
    LIMIT 10
  ) t;
  RETURN result;
END;
$$;

-- 12. Ingestion status
CREATE OR REPLACE FUNCTION public.pharmaco_ingestion_status()
RETURNS json LANGUAGE plpgsql SECURITY DEFINER AS $$
DECLARE result json;
BEGIN
  SELECT json_agg(row_to_json(t)) INTO result FROM (
    SELECT source, step, status, rows_inserted, started_at, completed_at
    FROM pharmaco_db.ingestion_log
    ORDER BY started_at
  ) t;
  RETURN result;
END;
$$;
