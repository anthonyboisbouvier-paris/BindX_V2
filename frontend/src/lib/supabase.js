import { createClient } from '@supabase/supabase-js'

// V9: read from env vars, fallback to hardcoded values for PharmacoDB compat
const SUPABASE_URL = import.meta.env.VITE_SUPABASE_URL || 'https://webkntghfzscrnuixfba.supabase.co'
const SUPABASE_KEY = import.meta.env.VITE_SUPABASE_ANON_KEY || 'sb_publishable_Zy18OaVJ4k_DgaafyiYaZA_9jB0ATmY'

export const supabase = createClient(SUPABASE_URL, SUPABASE_KEY)
