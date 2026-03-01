-- BindX V9 — Run logs table
-- Stores per-run log entries for frontend terminal display.

CREATE TABLE IF NOT EXISTS run_logs (
    id          UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    run_id      UUID NOT NULL REFERENCES runs(id) ON DELETE CASCADE,
    level       TEXT NOT NULL DEFAULT 'info',   -- info | warn | error | success
    message     TEXT NOT NULL,
    created_at  TIMESTAMPTZ NOT NULL DEFAULT now()
);

CREATE INDEX IF NOT EXISTS idx_run_logs_run_id ON run_logs(run_id);
CREATE INDEX IF NOT EXISTS idx_run_logs_created_at ON run_logs(run_id, created_at);
