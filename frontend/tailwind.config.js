/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // ═══ BindX Design System V9 ═══
        // Dark surfaces
        'bx-bg':       '#04060b',
        'bx-surface':  '#0f131d',
        'bx-elevated': '#161b2a',
        'bx-s1':       '#0a0d15',  // legacy alias
        'bx-s2':       '#0f131d',  // legacy alias
        'bx-s3':       '#141925',  // legacy alias
        // Dark text
        'bx-text':     '#eaedf4',
        'bx-sub':      '#8a93a8',
        'bx-dim':      '#5d6579',
        // Brand
        'bx-mint':     '#00e6a0',
        'bx-mint-dim': '#00c78a',
        'bx-cyan':     '#06b6d4',
        'bx-blue':     '#3b82f6',
        'bx-violet':   '#8b5cf6',
        // Semantic
        'bx-amber':    '#f59e0b',
        'bx-red':      '#ef4444',
        'bx-rose':     '#f43f5e',
        'bx-success':  '#10b981',
        // Light surfaces
        'bx-light-bg':      '#f4f5f8',
        'bx-light':         '#ffffff',
        'bx-light-warm':    '#f9fafb',
        'bx-light-border':  '#e2e5ec',
        'bx-light-border-s':'#edf0f4',
        // Light text
        'bx-light-text':    '#1e293b',
        'bx-light-text2':   '#475569',
        'bx-light-muted':   '#94a3b8',
        // Legacy aliases
        'dockit-blue':       '#1e3a5f',
        'dockit-blue-light': '#2d5a8e',
        'dockit-blue-dark':  '#132840',
        'dockit-green':      '#00e6a0',
        'dockit-green-dark': '#00c98b',
        'dockit-gray':       '#f4f5f8',
      },
      fontFamily: {
        sans:  ['Plus Jakarta Sans', 'Inter', 'system-ui', '-apple-system', 'sans-serif'],
        serif: ['Literata', 'Georgia', 'serif'],
        mono:  ['Fira Code', 'JetBrains Mono', 'monospace'],
      },
      borderRadius: {
        'card': '14px',
        'btn':  '9px',
        'badge': '5px',
        'input': '8px',
      },
      boxShadow: {
        'xs':         '0 1px 2px rgba(0,0,0,.03)',
        'card':       '0 1px 3px rgba(0,0,0,.04), 0 0 0 1px rgba(0,0,0,.02)',
        'card-hover': '0 4px 12px rgba(0,0,0,.05)',
        'md':         '0 4px 12px rgba(0,0,0,.05)',
        'lg':         '0 8px 24px rgba(0,0,0,.07)',
        'mint-glow':  '0 2px 8px rgba(0,230,160,.2)',
      },
      animation: {
        'pulse-slow': 'pulse 2s cubic-bezier(0.4, 0, 0.6, 1) infinite',
        'shimmer':    'shimmer 1.5s infinite',
        'rise':       'rise 0.5s ease-out both',
        'load-dot':   'ldp 1.8s ease-in-out infinite',
      },
      keyframes: {
        shimmer: {
          '0%':   { backgroundPosition: '-200% 0' },
          '100%': { backgroundPosition: '200% 0' },
        },
        rise: {
          from: { opacity: '0', transform: 'translateY(20px)' },
          to:   { opacity: '1', transform: 'translateY(0)' },
        },
        ldp: {
          '0%, 100%': { opacity: '.2', transform: 'scale(.7)' },
          '50%':      { opacity: '1', transform: 'scale(1)' },
        },
      },
    },
  },
  plugins: [],
}
