#!/bin/bash
# BindX V2 Development Shortcuts for Claude Code & Anthony
# Source this file: source tools/claude/shortcuts.sh

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

# =============================================================================
# NAVIGATION SHORTCUTS
# =============================================================================
alias bx="cd ~/BindX_V2"
alias bx-docs="cd ~/BindX_V2/docs && ls -la"
alias bx-backend="cd ~/BindX_V2/backend"
alias bx-frontend="cd ~/BindX_V2/frontend"
alias bx-tools="cd ~/BindX_V2/tools"
alias bx-scripts="cd ~/BindX_V2/scripts"

# =============================================================================
# DEVELOPMENT WORKFLOWS
# =============================================================================
alias bx-dev="cd ~/BindX_V2 && ./scripts/dev.sh"
alias bx-test="cd ~/BindX_V2 && ./scripts/test.sh"
alias bx-build="cd ~/BindX_V2 && ./scripts/build.sh"
alias bx-deploy="cd ~/BindX_V2 && ./scripts/deploy.sh"

# =============================================================================
# DOCKER MANAGEMENT
# =============================================================================
alias bx-up="cd ~/BindX_V2 && docker-compose up -d"
alias bx-down="cd ~/BindX_V2 && docker-compose down"
alias bx-restart="cd ~/BindX_V2 && docker-compose restart"
alias bx-rebuild="cd ~/BindX_V2 && docker-compose down && docker-compose up --build -d"
alias bx-logs="cd ~/BindX_V2 && docker-compose logs -f"
alias bx-logs-be="cd ~/BindX_V2 && docker-compose logs -f backend"
alias bx-logs-fe="cd ~/BindX_V2 && docker-compose logs -f frontend"

# Container shell access
alias bx-shell-be="cd ~/BindX_V2 && docker-compose exec backend bash"
alias bx-shell-fe="cd ~/BindX_V2 && docker-compose exec frontend sh"
alias bx-shell-db="cd ~/BindX_V2 && docker-compose exec postgres psql -U bindx -d bindx_dev"

# =============================================================================
# CLAUDE CODE INTEGRATION
# =============================================================================
alias claude-start="cd ~/BindX_V2 && claude"
alias claude-status="cd ~/BindX_V2 && git status && docker ps --format 'table {{.Names}}\t{{.Status}}\t{{.Ports}}'"
alias claude-health="cd ~/BindX_V2 && echo '🏥 Health Check:' && curl -f -s http://localhost:8000/health || echo '❌ Backend' && curl -f -s http://localhost:3000 >/dev/null && echo '✅ Frontend' || echo '❌ Frontend'"

# =============================================================================
# GIT SHORTCUTS
# =============================================================================
alias bx-status="cd ~/BindX_V2 && git status"
alias bx-diff="cd ~/BindX_V2 && git diff"
alias bx-add="cd ~/BindX_V2 && git add ."
alias bx-commit="cd ~/BindX_V2 && git add . && git commit -m"
alias bx-push="cd ~/BindX_V2 && git push"
alias bx-pull="cd ~/BindX_V2 && git pull"
alias bx-log="cd ~/BindX_V2 && git log --oneline -10"

# =============================================================================
# FILE MANAGEMENT
# =============================================================================
alias bx-clean="cd ~/BindX_V2 && docker system prune -f && docker volume prune -f"
alias bx-backup="cd ~/BindX_V2 && tar -czf ~/bindx-backup-$(date +%Y%m%d).tar.gz . --exclude=data --exclude=node_modules --exclude=.git"
alias bx-find="find ~/BindX_V2 -name"
alias bx-grep="cd ~/BindX_V2 && grep -r --exclude-dir=node_modules --exclude-dir=.git"

# =============================================================================
# DEVELOPMENT TOOLS
# =============================================================================
# Python environment
alias bx-py="cd ~/BindX_V2/backend && python3"
alias bx-pip="cd ~/BindX_V2/backend && pip3"
alias bx-pytest="cd ~/BindX_V2/backend && python3 -m pytest"

# Node.js environment  
alias bx-npm="cd ~/BindX_V2/frontend && npm"
alias bx-yarn="cd ~/BindX_V2/frontend && yarn"

# Database
alias bx-psql="psql -h localhost -U bindx -d bindx_dev"
alias bx-redis="redis-cli -h localhost"

# =============================================================================
# MOLECULAR TOOLS
# =============================================================================
alias bx-vina="docker run --rm -v ~/BindX_V2/data:/data autodocksuite/autodock-vina"
alias bx-obabel="docker run --rm -v ~/BindX_V2/data:/data openbabel/openbabel obabel"

# =============================================================================
# MONITORING
# =============================================================================
alias bx-ps="cd ~/BindX_V2 && docker-compose ps"
alias bx-top="cd ~/BindX_V2 && docker stats --format 'table {{.Container}}\t{{.CPUPerc}}\t{{.MemUsage}}'"
alias bx-ports="cd ~/BindX_V2 && docker-compose ps | grep -E '(frontend|backend|postgres|redis)'"

# =============================================================================
# HELPFUL FUNCTIONS
# =============================================================================

# Quick project status
bx-info() {
    echo -e "${BLUE}🧬 BindX V2 Project Information${NC}"
    echo "================================"
    echo -e "${YELLOW}📁 Location:${NC} $(pwd)"
    echo -e "${YELLOW}🌿 Branch:${NC} $(git branch --show-current 2>/dev/null || echo 'Not a git repo')"
    echo -e "${YELLOW}📊 Containers:${NC}"
    cd ~/BindX_V2 && docker-compose ps 2>/dev/null || echo "Docker not running"
    echo -e "${YELLOW}🌐 Services:${NC}"
    echo "   Frontend:  http://localhost:3000"
    echo "   Backend:   http://localhost:8000"
    echo "   API Docs:  http://localhost:8000/docs"
}

# Start complete development environment
bx-start() {
    echo -e "${GREEN}🚀 Starting BindX V2 development environment...${NC}"
    cd ~/BindX_V2
    ./scripts/dev.sh
}

# Stop everything
bx-stop() {
    echo -e "${YELLOW}🛑 Stopping BindX V2 services...${NC}"
    cd ~/BindX_V2
    docker-compose down
}

# Reset development environment
bx-reset() {
    echo -e "${YELLOW}🔄 Resetting BindX V2 development environment...${NC}"
    cd ~/BindX_V2
    docker-compose down
    docker-compose up --build -d
}

# =============================================================================
# CLAUDE CODE HELPERS  
# =============================================================================

# Load common Claude tasks
claude-tasks() {
    cat ~/BindX_V2/tools/claude/tasks.md
}

# Quick Claude status
claude-ready() {
    cd ~/BindX_V2
    echo -e "${BLUE}🤖 Claude Code Environment Check${NC}"
    echo "================================="
    echo -e "${YELLOW}📁 Project:${NC} $(pwd)"
    echo -e "${YELLOW}🌿 Branch:${NC} $(git branch --show-current)"
    echo -e "${YELLOW}📊 Containers:${NC} $(docker-compose ps | grep -c 'Up')/$(docker-compose ps | wc -l) running"
    echo -e "${YELLOW}🔧 Ready:${NC} $([ -f .env ] && echo '✅ .env configured' || echo '❌ Need .env setup')"
}

# =============================================================================
# LOAD CONFIRMATION
# =============================================================================

echo -e "${GREEN}🤖 BindX V2 shortcuts loaded successfully!${NC}"
echo ""
echo -e "${BLUE}📁 Navigation:${NC}"
echo "   bx           → go to project root"
echo "   bx-docs      → go to documentation"  
echo "   bx-backend   → go to backend code"
echo "   bx-frontend  → go to frontend code"
echo ""
echo -e "${BLUE}🚀 Development:${NC}"
echo "   bx-dev       → start development environment"
echo "   bx-test      → run all tests"
echo "   bx-build     → build all containers"
echo "   bx-logs      → view all container logs"
echo ""
echo -e "${BLUE}🐳 Docker:${NC}"
echo "   bx-up        → start containers"
echo "   bx-down      → stop containers"
echo "   bx-restart   → restart containers"
echo "   bx-rebuild   → rebuild & restart containers"
echo ""
echo -e "${BLUE}🤖 Claude Code:${NC}"
echo "   claude-start → start AI coding session"
echo "   claude-status → project & container status"
echo "   claude-health → health check all services"
echo ""
echo -e "${BLUE}ℹ️  Information:${NC}"
echo "   bx-info      → project overview"
echo "   bx-start     → complete startup"
echo "   bx-reset     → reset environment"
echo ""
echo -e "${YELLOW}💡 Usage: Run 'bx-start' to begin development${NC}"
