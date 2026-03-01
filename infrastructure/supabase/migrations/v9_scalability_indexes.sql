-- BindX V9 — Scalability indexes for 100K+ molecules
-- Run against Supabase PostgreSQL

-- Composite index for paginated listing (keyset pagination on created_at)
CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_molecules_phase_created
  ON molecules(phase_id, created_at DESC);

-- Partial index for bookmarked-only queries
CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_molecules_phase_bookmarked
  ON molecules(phase_id, bookmarked) WHERE bookmarked = true;

-- Index for name search/sort
CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_molecules_phase_name
  ON molecules(phase_id, name);

-- Index for property lookups (lateral join for sort-by-property)
CREATE INDEX CONCURRENTLY IF NOT EXISTS idx_mol_props_mol_name
  ON molecule_properties(molecule_id, property_name);
