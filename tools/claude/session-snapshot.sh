#!/bin/bash
# Session snapshot system for Claude Code

SNAPSHOT_DIR="$HOME/.claude-snapshots/$(date +%Y%m%d)"

case "$1" in
    "start")
        mkdir -p "$SNAPSHOT_DIR"
        echo "📸 Creating pre-session snapshot..."
        
        # Code metrics
        find . -name "*.py" -o -name "*.js" -o -name "*.jsx" | xargs wc -l > "$SNAPSHOT_DIR/lines-before.txt" 2>/dev/null
        git log --oneline -10 > "$SNAPSHOT_DIR/git-before.txt" 2>/dev/null
        docker ps > "$SNAPSHOT_DIR/containers-before.txt" 2>/dev/null
        
        echo "✅ Pre-session snapshot saved"
        ;;
        
    "end")
        echo "📸 Creating post-session snapshot..."
        
        # Code metrics
        find . -name "*.py" -o -name "*.js" -o -name "*.jsx" | xargs wc -l > "$SNAPSHOT_DIR/lines-after.txt" 2>/dev/null
        git log --oneline -10 > "$SNAPSHOT_DIR/git-after.txt" 2>/dev/null
        
        # Session summary
        echo "📊 Session Summary:" > "$SNAPSHOT_DIR/summary.txt"
        echo "==================" >> "$SNAPSHOT_DIR/summary.txt"
        echo "Date: $(date)" >> "$SNAPSHOT_DIR/summary.txt"
        echo "Project: $(pwd)" >> "$SNAPSHOT_DIR/summary.txt"
        echo "Branch: $(git branch --show-current 2>/dev/null || echo 'no-git')" >> "$SNAPSHOT_DIR/summary.txt"
        
        # Git commits during session
        commits_today=$(git log --oneline --since="1 day ago" | wc -l)
        echo "Commits today: $commits_today" >> "$SNAPSHOT_DIR/summary.txt"
        
        cat "$SNAPSHOT_DIR/summary.txt"
        echo "✅ Post-session snapshot saved"
        ;;
        
    *)
        echo "Usage: session-snapshot.sh {start|end}"
        ;;
esac
