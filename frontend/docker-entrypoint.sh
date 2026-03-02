#!/bin/sh
# Custom entrypoint: switch nginx config based on environment
# - BACKEND_HOST set → production mode (Render) with basic auth + proxy to cloud backend
# - BACKEND_HOST unset → local mode (docker-compose) with proxy to http://backend:8000

if [ -n "$BACKEND_HOST" ]; then
  export BACKEND_URL="https://${BACKEND_HOST}"
  echo "Production mode: proxying /api/ to $BACKEND_URL"
  envsubst '${BACKEND_URL} ${PORT}' < /etc/nginx/nginx.prod.conf.template > /etc/nginx/conf.d/default.conf
else
  echo "Local mode: using docker-compose config (proxy to http://backend:8000)"
fi

exec "$@"
