FROM node:20-alpine AS build
WORKDIR /app

# Install dependencies first (layer caching)
COPY package*.json ./
RUN npm ci || npm install

# Copy source and build
COPY . .
RUN npm run build

# Production image with nginx
FROM nginx:alpine
RUN apk add --no-cache curl apache2-utils

# Generate htpasswd for basic auth (user: bindx, password: test1234)
RUN htpasswd -cb /etc/nginx/.htpasswd bindx test1234

# Local config (docker-compose) — default
COPY nginx.conf /etc/nginx/conf.d/default.conf

# Production template (Render) — activated by BACKEND_URL env var
COPY nginx.prod.conf.template /etc/nginx/nginx.prod.conf.template

# Custom entrypoint to switch config based on environment
COPY docker-entrypoint.sh /docker-entrypoint.sh
RUN chmod +x /docker-entrypoint.sh

COPY --from=build /app/dist /usr/share/nginx/html

# Remove default nginx entrypoint template processing to avoid conflicts
RUN rm -rf /etc/nginx/templates

HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
    CMD curl -sf http://localhost:${PORT:-80}/ || exit 1

EXPOSE ${PORT:-80}
ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["nginx", "-g", "daemon off;"]
