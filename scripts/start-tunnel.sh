#!/bin/bash
# BindX — Start backend + tunnel
# Run this once when you boot your PC. The site stays online as long as this runs.
# URL: https://anthonyboisbouvier-paris.github.io/BindX_V2/
# API:  https://bindx-api.loca.lt

set -e

echo "=== BindX — Starting backend + tunnel ==="

# 1. Start Docker containers (if not already running)
if ! docker ps | grep -q bindx_v2-backend; then
  echo "[1/2] Starting Docker containers..."
  cd ~/BindX_V2
  docker-compose up -d
  echo "  Waiting for backend to be healthy..."
  sleep 10
else
  echo "[1/2] Docker containers already running"
fi

# 2. Start localtunnel with fixed subdomain
echo "[2/2] Starting tunnel (bindx-api.loca.lt)..."
pkill -f localtunnel 2>/dev/null || true
sleep 1

npx localtunnel --port 8000 --subdomain bindx-api &
LT_PID=$!
sleep 5

# Verify
if curl -s https://bindx-api.loca.lt/api/v9/health -H "bypass-tunnel-reminder: true" | grep -q '"ok"'; then
  echo ""
  echo "=== ONLINE ==="
  echo "Site:    https://anthonyboisbouvier-paris.github.io/BindX_V2/"
  echo "API:     https://bindx-api.loca.lt"
  echo "Login:   anthony.boisbouvier@gmail.com / test1234"
  echo ""
  echo "Press Ctrl+C to stop the tunnel."
  wait $LT_PID
else
  echo "ERROR: API not reachable via tunnel"
  kill $LT_PID 2>/dev/null
  exit 1
fi
