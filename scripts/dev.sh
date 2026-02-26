#!/bin/bash
# BindX V2 Development Setup Script
# Run this to start the complete development environment

set -e

echo "🧬 BindX V2 Development Environment Setup"
echo "========================================"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Check if we're in the right directory
if [ ! -f "README.md" ] || ! grep -q "BindX V2" README.md 2>/dev/null; then
    echo -e "${RED}❌ Error: Run this script from the BindX_V2 root directory${NC}"
    echo "Current directory: $(pwd)"
    exit 1
fi

echo -e "${BLUE}📁 Current directory: $(pwd)${NC}"

# Check prerequisites
echo -e "\n${YELLOW}🔍 Checking prerequisites...${NC}"

# Docker
if command -v docker >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Docker found: $(docker --version | cut -d' ' -f3)${NC}"
else
    echo -e "${RED}❌ Docker not found. Please install Docker Desktop${NC}"
    echo "Download: https://www.docker.com/products/docker-desktop"
    exit 1
fi

# Docker Compose
if command -v docker-compose >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Docker Compose found: $(docker-compose --version | cut -d' ' -f3)${NC}"
else
    echo -e "${RED}❌ Docker Compose not found${NC}"
    exit 1
fi

# Node.js (optional check)
if command -v node >/dev/null 2>&1; then
    NODE_VERSION=$(node --version)
    echo -e "${GREEN}✅ Node.js found: $NODE_VERSION${NC}"
    # Check if version is 18+
    if [ "${NODE_VERSION#v}" \< "18" ]; then
        echo -e "${YELLOW}⚠️  Node.js 18+ recommended (current: $NODE_VERSION)${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  Node.js not found in PATH (not required for Docker development)${NC}"
fi

# Python (optional check)
if command -v python3 >/dev/null 2>&1; then
    PYTHON_VERSION=$(python3 --version 2>&1 | cut -d' ' -f2)
    echo -e "${GREEN}✅ Python found: $PYTHON_VERSION${NC}"
else
    echo -e "${YELLOW}⚠️  Python3 not found in PATH (not required for Docker development)${NC}"
fi

# Environment configuration
echo -e "\n${YELLOW}📝 Setting up environment...${NC}"

if [ ! -f ".env" ]; then
    echo -e "${BLUE}Creating .env from template...${NC}"
    cp .env.example .env
    echo -e "${YELLOW}⚠️  IMPORTANT: Edit .env with your Supabase credentials!${NC}"
    echo -e "${BLUE}ℹ️  See docs/DEVELOPMENT.md for detailed setup guide${NC}"
else
    echo -e "${GREEN}✅ .env file already exists${NC}"
fi

# Create data directories (gitignored)
echo -e "\n${YELLOW}📂 Creating data directories...${NC}"
mkdir -p data/{structures,libraries,models,cache,artifacts}
mkdir -p logs
echo -e "${GREEN}✅ Data directories created${NC}"

# Docker services
echo -e "\n${YELLOW}🐳 Starting Docker services...${NC}"

# Check if docker-compose.yml exists
if [ ! -f "docker-compose.yml" ]; then
    echo -e "${YELLOW}⚠️  docker-compose.yml not found, creating basic configuration...${NC}"
    
    cat > docker-compose.yml << 'EOF'
version: '3.8'

services:
  # Redis for Celery
  redis:
    image: redis:7-alpine
    ports:
      - "6379:6379"
    volumes:
      - redis_data:/data
    restart: unless-stopped

  # PostgreSQL (local development - use Supabase in production)
  postgres:
    image: postgres:15
    environment:
      POSTGRES_DB: bindx_dev
      POSTGRES_USER: bindx
      POSTGRES_PASSWORD: bindx_dev_password
    ports:
      - "5432:5432"
    volumes:
      - postgres_data:/var/lib/postgresql/data
    restart: unless-stopped

  # Backend API (FastAPI)
  backend:
    build: 
      context: ./backend
      dockerfile: Dockerfile
    ports:
      - "8000:8000"
    environment:
      - DATABASE_URL=postgresql://bindx:bindx_dev_password@postgres:5432/bindx_dev
      - REDIS_URL=redis://redis:6379/0
    volumes:
      - ./backend:/app
      - ./data:/data
    depends_on:
      - postgres
      - redis
    restart: unless-stopped

  # Frontend (React)
  frontend:
    build:
      context: ./frontend
      dockerfile: Dockerfile
    ports:
      - "3000:3000"
    volumes:
      - ./frontend:/app
      - /app/node_modules
    environment:
      - REACT_APP_API_URL=http://localhost:8000
    depends_on:
      - backend
    restart: unless-stopped

  # Celery Worker
  celery:
    build:
      context: ./backend
      dockerfile: Dockerfile
    command: celery -A app.celery_app worker --loglevel=info
    environment:
      - DATABASE_URL=postgresql://bindx:bindx_dev_password@postgres:5432/bindx_dev
      - REDIS_URL=redis://redis:6379/0
    volumes:
      - ./backend:/app
      - ./data:/data
    depends_on:
      - postgres
      - redis
    restart: unless-stopped

volumes:
  postgres_data:
  redis_data:
EOF
    echo -e "${GREEN}✅ Basic docker-compose.yml created${NC}"
fi

# Start services
echo -e "${BLUE}Starting Docker Compose services...${NC}"
docker-compose up -d

# Wait for services to be ready
echo -e "\n${YELLOW}⏳ Waiting for services to initialize...${NC}"
sleep 10

# Health checks
echo -e "\n${YELLOW}🏥 Performing health checks...${NC}"

# Check backend
if curl -f -s http://localhost:8000/health >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Backend API healthy${NC}"
elif curl -f -s http://localhost:8000 >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Backend API responding${NC}"
else
    echo -e "${YELLOW}⚠️  Backend API not ready yet (this is normal during first startup)${NC}"
fi

# Check frontend
if curl -f -s http://localhost:3000 >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Frontend responding${NC}"
else
    echo -e "${YELLOW}⚠️  Frontend not ready yet (this is normal during first startup)${NC}"
fi

# Check Redis
if redis-cli -h localhost ping >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Redis responding${NC}"
elif nc -z localhost 6379 >/dev/null 2>&1; then
    echo -e "${GREEN}✅ Redis port accessible${NC}"
else
    echo -e "${YELLOW}⚠️  Redis not accessible${NC}"
fi

# Final status
echo -e "\n${GREEN}🎉 BindX V2 Development Environment Ready!${NC}"
echo -e "========================================"
echo -e "${BLUE}🌐 Frontend:${NC}  http://localhost:3000"
echo -e "${BLUE}🐍 Backend:${NC}   http://localhost:8000"
echo -e "${BLUE}📚 API Docs:${NC}  http://localhost:8000/docs"
echo -e "${BLUE}📊 Redis:${NC}     localhost:6379"
echo -e "${BLUE}🗄️  Postgres:${NC} localhost:5432"

# Development tips
echo -e "\n${YELLOW}💡 Development Tips:${NC}"
echo -e "${BLUE}View logs:${NC}     docker-compose logs -f"
echo -e "${BLUE}Restart:${NC}      docker-compose restart"
echo -e "${BLUE}Stop all:${NC}     docker-compose down"
echo -e "${BLUE}Clean rebuild:${NC} docker-compose down && docker-compose up --build -d"

# Claude Code integration
echo -e "\n${YELLOW}🤖 Claude Code Integration:${NC}"
echo -e "${BLUE}Start AI coding:${NC} claude"
echo -e "${BLUE}Load shortcuts:${NC} source tools/claude/shortcuts.sh"

# Next steps
echo -e "\n${YELLOW}📋 Next Steps:${NC}"
echo "1. ${BLUE}Configure Supabase${NC}: Edit .env with your Supabase credentials"
echo "2. ${BLUE}Review docs${NC}: Check docs/DEVELOPMENT.md for detailed setup"
echo "3. ${BLUE}Start coding${NC}: Open http://localhost:3000 and begin development"
echo "4. ${BLUE}Use Claude Code${NC}: Run 'claude' for AI-assisted development"

echo -e "\n${GREEN}Happy coding! 🚀${NC}"
