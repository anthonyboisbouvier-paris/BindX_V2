-- ============================================================================
-- BindX V9 — Database Schema
-- Target: Supabase Cloud (PostgreSQL 15+)
-- Schema: public (no conflict with pharmaco_db.*)
--
-- Run this in Supabase SQL Editor.
-- ============================================================================

-- 1. projects
CREATE TABLE projects (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID REFERENCES auth.users(id),
  name TEXT NOT NULL,
  description TEXT,
  -- Target info
  target_input_type TEXT CHECK (target_input_type IN ('uniprot', 'sequence', 'pdb')),
  target_input_value TEXT NOT NULL,
  target_name TEXT,
  target_pdb_id TEXT,
  target_preview JSONB,
  -- Structure source info
  structure_source TEXT,
  structure_resolution FLOAT,
  structure_method TEXT,
  cocrystal_ligand TEXT,
  -- Pockets detected
  pockets_detected JSONB,
  -- ChEMBL info
  chembl_actives_count INTEGER,
  chembl_median_ic50 FLOAT,
  --
  status TEXT DEFAULT 'active',
  notification_email TEXT,
  created_at TIMESTAMPTZ DEFAULT now(),
  updated_at TIMESTAMPTZ DEFAULT now()
);

-- 2. campaigns
CREATE TABLE campaigns (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  project_id UUID REFERENCES projects(id) ON DELETE CASCADE,
  name TEXT NOT NULL,
  pocket_config JSONB,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- 3. phases
CREATE TABLE phases (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  campaign_id UUID REFERENCES campaigns(id) ON DELETE CASCADE,
  type TEXT NOT NULL CHECK (type IN ('hit_discovery', 'hit_to_lead', 'lead_optimization')),
  status TEXT DEFAULT 'active' CHECK (status IN ('active', 'frozen', 'completed')),
  frozen_at TIMESTAMPTZ,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- 4. runs
CREATE TABLE runs (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  phase_id UUID REFERENCES phases(id) ON DELETE CASCADE,
  type TEXT NOT NULL CHECK (type IN ('import', 'calculation', 'generation')),
  calculation_types TEXT[],
  status TEXT DEFAULT 'created' CHECK (status IN ('created', 'queued', 'running', 'completed', 'failed')),
  config JSONB NOT NULL,
  input_molecule_ids UUID[],
  input_source TEXT,
  input_file_path TEXT,
  progress INTEGER DEFAULT 0,
  current_step TEXT,
  estimated_duration_seconds INTEGER,
  started_at TIMESTAMPTZ,
  completed_at TIMESTAMPTZ,
  error_message TEXT,
  archived BOOLEAN DEFAULT false,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- 5. molecules
CREATE TABLE molecules (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  phase_id UUID REFERENCES phases(id) ON DELETE CASCADE,
  smiles TEXT NOT NULL,
  canonical_smiles TEXT NOT NULL,
  name TEXT,
  source_run_id UUID REFERENCES runs(id),
  bookmarked BOOLEAN DEFAULT false,
  generation_level INTEGER DEFAULT 0,
  parent_molecule_id UUID REFERENCES molecules(id),
  ai_generated BOOLEAN DEFAULT false,
  created_at TIMESTAMPTZ DEFAULT now(),
  UNIQUE(phase_id, canonical_smiles)
);

-- 6. molecule_properties
CREATE TABLE molecule_properties (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  molecule_id UUID REFERENCES molecules(id) ON DELETE CASCADE,
  run_id UUID REFERENCES runs(id),
  property_name TEXT NOT NULL,
  property_value JSONB NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now(),
  UNIQUE(molecule_id, property_name, run_id)
);

-- 7. calculation_cache
CREATE TABLE calculation_cache (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  cache_key TEXT UNIQUE NOT NULL,
  result JSONB NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- 8. artifacts
CREATE TABLE artifacts (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  run_id UUID REFERENCES runs(id) ON DELETE CASCADE,
  type TEXT NOT NULL,
  storage_path TEXT NOT NULL,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- 9. audit_log
CREATE TABLE audit_log (
  id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
  user_id UUID REFERENCES auth.users(id),
  action TEXT NOT NULL,
  entity_type TEXT,
  entity_id UUID,
  details JSONB,
  created_at TIMESTAMPTZ DEFAULT now()
);

-- ============================================================================
-- Indexes
-- ============================================================================

CREATE INDEX idx_projects_user ON projects(user_id);
CREATE INDEX idx_campaigns_project ON campaigns(project_id);
CREATE INDEX idx_phases_campaign ON phases(campaign_id);
CREATE INDEX idx_runs_phase ON runs(phase_id);
CREATE INDEX idx_molecules_phase ON molecules(phase_id);
CREATE INDEX idx_molecules_canonical ON molecules(phase_id, canonical_smiles);
CREATE INDEX idx_molecules_parent ON molecules(parent_molecule_id);
CREATE INDEX idx_mol_props_molecule ON molecule_properties(molecule_id);
CREATE INDEX idx_mol_props_name ON molecule_properties(property_name);
CREATE INDEX idx_cache_key ON calculation_cache(cache_key);
CREATE INDEX idx_audit_created ON audit_log(created_at);

-- ============================================================================
-- Trigger: auto-update updated_at on projects
-- ============================================================================

CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
  NEW.updated_at = now();
  RETURN NEW;
END;
$$ LANGUAGE plpgsql;

CREATE TRIGGER trg_projects_updated_at
  BEFORE UPDATE ON projects
  FOR EACH ROW
  EXECUTE FUNCTION update_updated_at_column();
