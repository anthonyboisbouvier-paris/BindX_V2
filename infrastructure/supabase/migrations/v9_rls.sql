-- ============================================================================
-- BindX V9 — Row Level Security Policies
-- Target: Supabase Cloud (PostgreSQL 15+)
--
-- Run this in Supabase SQL Editor AFTER v9_schema.sql.
-- ============================================================================

-- ============================================================================
-- Enable RLS on all V9 tables
-- ============================================================================

ALTER TABLE projects ENABLE ROW LEVEL SECURITY;
ALTER TABLE campaigns ENABLE ROW LEVEL SECURITY;
ALTER TABLE phases ENABLE ROW LEVEL SECURITY;
ALTER TABLE runs ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecules ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecule_properties ENABLE ROW LEVEL SECURITY;
ALTER TABLE calculation_cache ENABLE ROW LEVEL SECURITY;
ALTER TABLE artifacts ENABLE ROW LEVEL SECURITY;
ALTER TABLE audit_log ENABLE ROW LEVEL SECURITY;

-- ============================================================================
-- projects: direct ownership
-- ============================================================================

CREATE POLICY "projects_select_own" ON projects
  FOR SELECT USING (auth.uid() = user_id);

CREATE POLICY "projects_insert_own" ON projects
  FOR INSERT WITH CHECK (auth.uid() = user_id);

CREATE POLICY "projects_update_own" ON projects
  FOR UPDATE USING (auth.uid() = user_id);

CREATE POLICY "projects_delete_own" ON projects
  FOR DELETE USING (auth.uid() = user_id);

-- ============================================================================
-- campaigns: transitive ownership via projects
-- ============================================================================

CREATE POLICY "campaigns_select_own" ON campaigns
  FOR SELECT USING (
    EXISTS (SELECT 1 FROM projects WHERE projects.id = campaigns.project_id AND projects.user_id = auth.uid())
  );

CREATE POLICY "campaigns_insert_own" ON campaigns
  FOR INSERT WITH CHECK (
    EXISTS (SELECT 1 FROM projects WHERE projects.id = campaigns.project_id AND projects.user_id = auth.uid())
  );

CREATE POLICY "campaigns_update_own" ON campaigns
  FOR UPDATE USING (
    EXISTS (SELECT 1 FROM projects WHERE projects.id = campaigns.project_id AND projects.user_id = auth.uid())
  );

CREATE POLICY "campaigns_delete_own" ON campaigns
  FOR DELETE USING (
    EXISTS (SELECT 1 FROM projects WHERE projects.id = campaigns.project_id AND projects.user_id = auth.uid())
  );

-- ============================================================================
-- phases: transitive via campaigns → projects
-- ============================================================================

CREATE POLICY "phases_select_own" ON phases
  FOR SELECT USING (
    EXISTS (
      SELECT 1 FROM campaigns
      JOIN projects ON projects.id = campaigns.project_id
      WHERE campaigns.id = phases.campaign_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "phases_insert_own" ON phases
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM campaigns
      JOIN projects ON projects.id = campaigns.project_id
      WHERE campaigns.id = phases.campaign_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "phases_update_own" ON phases
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM campaigns
      JOIN projects ON projects.id = campaigns.project_id
      WHERE campaigns.id = phases.campaign_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "phases_delete_own" ON phases
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM campaigns
      JOIN projects ON projects.id = campaigns.project_id
      WHERE campaigns.id = phases.campaign_id AND projects.user_id = auth.uid()
    )
  );

-- ============================================================================
-- runs: transitive via phases → campaigns → projects
-- ============================================================================

CREATE POLICY "runs_select_own" ON runs
  FOR SELECT USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = runs.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "runs_insert_own" ON runs
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = runs.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "runs_update_own" ON runs
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = runs.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "runs_delete_own" ON runs
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = runs.phase_id AND projects.user_id = auth.uid()
    )
  );

-- ============================================================================
-- molecules: transitive via phases → campaigns → projects
-- ============================================================================

CREATE POLICY "molecules_select_own" ON molecules
  FOR SELECT USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = molecules.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "molecules_insert_own" ON molecules
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = molecules.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "molecules_update_own" ON molecules
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = molecules.phase_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "molecules_delete_own" ON molecules
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM phases
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE phases.id = molecules.phase_id AND projects.user_id = auth.uid()
    )
  );

-- ============================================================================
-- molecule_properties: transitive via molecules → phases → campaigns → projects
-- ============================================================================

CREATE POLICY "mol_props_select_own" ON molecule_properties
  FOR SELECT USING (
    EXISTS (
      SELECT 1 FROM molecules
      JOIN phases ON phases.id = molecules.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE molecules.id = molecule_properties.molecule_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "mol_props_insert_own" ON molecule_properties
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM molecules
      JOIN phases ON phases.id = molecules.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE molecules.id = molecule_properties.molecule_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "mol_props_update_own" ON molecule_properties
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM molecules
      JOIN phases ON phases.id = molecules.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE molecules.id = molecule_properties.molecule_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "mol_props_delete_own" ON molecule_properties
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM molecules
      JOIN phases ON phases.id = molecules.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE molecules.id = molecule_properties.molecule_id AND projects.user_id = auth.uid()
    )
  );

-- ============================================================================
-- artifacts: transitive via runs → phases → campaigns → projects
-- ============================================================================

CREATE POLICY "artifacts_select_own" ON artifacts
  FOR SELECT USING (
    EXISTS (
      SELECT 1 FROM runs
      JOIN phases ON phases.id = runs.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE runs.id = artifacts.run_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "artifacts_insert_own" ON artifacts
  FOR INSERT WITH CHECK (
    EXISTS (
      SELECT 1 FROM runs
      JOIN phases ON phases.id = runs.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE runs.id = artifacts.run_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "artifacts_update_own" ON artifacts
  FOR UPDATE USING (
    EXISTS (
      SELECT 1 FROM runs
      JOIN phases ON phases.id = runs.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE runs.id = artifacts.run_id AND projects.user_id = auth.uid()
    )
  );

CREATE POLICY "artifacts_delete_own" ON artifacts
  FOR DELETE USING (
    EXISTS (
      SELECT 1 FROM runs
      JOIN phases ON phases.id = runs.phase_id
      JOIN campaigns ON campaigns.id = phases.campaign_id
      JOIN projects ON projects.id = campaigns.project_id
      WHERE runs.id = artifacts.run_id AND projects.user_id = auth.uid()
    )
  );

-- ============================================================================
-- audit_log: direct ownership
-- ============================================================================

CREATE POLICY "audit_log_select_own" ON audit_log
  FOR SELECT USING (auth.uid() = user_id);

CREATE POLICY "audit_log_insert_own" ON audit_log
  FOR INSERT WITH CHECK (auth.uid() = user_id);

-- ============================================================================
-- calculation_cache: read for authenticated, write for service_role only
-- ============================================================================

CREATE POLICY "cache_select_authenticated" ON calculation_cache
  FOR SELECT USING (auth.role() = 'authenticated');

-- No INSERT/UPDATE/DELETE policies for authenticated users.
-- Only service_role (backend) can write to calculation_cache (bypasses RLS).
