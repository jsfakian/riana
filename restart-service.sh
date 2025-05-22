#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------
# start-django.sh — activate venv, stop & restart
# -----------------------------------------------

# Go to project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
cd "$SCRIPT_DIR"

# Activate it
# shellcheck source=/dev/null
source venv/bin/activate
source .env

# Stop any running Django dev server on port 8000
PORT=8001

# Find PIDs listening on $PORT
PIDS="$(lsof -ti tcp:"$PORT" || true)"
if [ -n "$PIDS" ]; then
  echo "→ Stopping existing process on port $PORT (PIDs: $PIDS)"
  kill $PIDS
  # wait for them to actually exit
  for pid in $PIDS; do
    while kill -0 "$pid" &>/dev/null; do
      sleep 0.1
    done
  done
fi

# Launch Django dev server
echo "→ Starting Django dev server on 0.0.0.0:$PORT"
exec nohup python manage.py runserver 0.0.0.0:"$PORT" &

