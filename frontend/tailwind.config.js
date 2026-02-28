/** @type {import('tailwindcss').Config} */
export default {
  content: [
    "./index.html",
    "./src/**/*.{js,ts,jsx,tsx}",
  ],
  theme: {
    extend: {
      colors: {
        // BindX design system (v10)
        'bx-bg':     '#04060b',
        'bx-s1':     '#0a0d15',
        'bx-s2':     '#0f131d',
        'bx-s3':     '#141925',
        'bx-mint':   '#00e6a0',
        'bx-blue':   '#3b82f6',
        'bx-violet': '#8b5cf6',
        'bx-amber':  '#f59e0b',
        'bx-rose':   '#f43f5e',
        'bx-cyan':   '#06b6d4',
        'bx-text':   '#eaedf4',
        'bx-sub':    '#8a93a8',
        'bx-dim':    '#5d6579',
        // Legacy aliases (for existing app components)
        'dockit-blue':       '#1e3a5f',
        'dockit-blue-light': '#2d5a8e',
        'dockit-blue-dark':  '#132840',
        'dockit-green':      '#00e6a0',
        'dockit-green-dark': '#00c98b',
        'dockit-gray':       '#f8fafc',
      },
      fontFamily: {
        sans:  ['Plus Jakarta Sans', 'Inter', 'system-ui', '-apple-system', 'sans-serif'],
        serif: ['Literata', 'Georgia', 'serif'],
        mono:  ['Fira Code', 'JetBrains Mono', 'monospace'],
      },
      boxShadow: {
        'card':       '0 4px 6px -1px rgba(0,0,0,0.1), 0 2px 4px -1px rgba(0,0,0,0.06)',
        'card-hover': '0 10px 15px -3px rgba(0,0,0,0.1), 0 4px 6px -2px rgba(0,0,0,0.05)',
      },
      animation: {
        'pulse-slow': 'pulse 2s cubic-bezier(0.4, 0, 0.6, 1) infinite',
        'shimmer':    'shimmer 1.5s infinite',
        'rise':       'rise 0.5s ease-out both',
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
      },
    },
  },
  plugins: [],
}
