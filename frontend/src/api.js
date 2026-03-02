import axios from 'axios'

const BASE_URL = import.meta.env.VITE_API_URL || '/api'

const apiClient = axios.create({
  baseURL: BASE_URL,
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
    // Bypass ngrok free-tier browser warning on API calls
    ...(BASE_URL.includes('ngrok') ? { 'ngrok-skip-browser-warning': 'true' } : {}),
  },
})

// Request interceptor for logging
apiClient.interceptors.request.use(
  (config) => {
    console.debug(`[API] ${config.method?.toUpperCase()} ${config.url}`)
    return config
  },
  (error) => Promise.reject(error)
)

// Response interceptor for error handling
apiClient.interceptors.response.use(
  (response) => response,
  (error) => {
    if (error.code === 'ECONNREFUSED' || error.code === 'ERR_NETWORK') {
      console.error('[API] Backend unreachable')
      error.userMessage = 'The backend server is unreachable. Make sure the backend is running on port 8000.'
    } else if (error.response) {
      const status = error.response.status
      if (status === 404) {
        error.userMessage = 'Resource not found.'
      } else if (status === 422) {
        const detail = error.response.data?.detail
        if (Array.isArray(detail)) {
          error.userMessage = detail.map(d => d.msg).join(', ')
        } else {
          error.userMessage = detail || 'Invalid data.'
        }
      } else if (status >= 500) {
        error.userMessage = 'Internal server error. Check the backend logs.'
      } else {
        error.userMessage = error.response.data?.detail || 'An error occurred.'
      }
    } else {
      error.userMessage = 'Connection error with the server.'
    }
    return Promise.reject(error)
  }
)

// ---------------------------------------------------------------------------
// Legacy V8 endpoints (still used by V9 frontend)
// ---------------------------------------------------------------------------

export function getReportUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/report`
}

export function getDownloadUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/download`
}

export function getProteinUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/protein`
}

export function getPoseUrl(jobId, index) {
  return `${BASE_URL}/jobs/${jobId}/pose/${index}`
}

export async function previewTarget(uniprotId) {
  const response = await apiClient.post('/preview-target', { uniprot_id: uniprotId }, { timeout: 120000 })
  return response.data
}

export async function detectPockets(downloadUrl, label, ligandId) {
  const response = await apiClient.post('/detect-pockets', {
    download_url: downloadUrl,
    label: label || 'structure',
    ligand_id: ligandId || null,
  }, { timeout: 120000 })
  return response.data
}

export async function analyzeScaffold(smiles) {
  const response = await apiClient.post('/molecule/analyze-scaffold', { smiles })
  return response.data
}

export async function queryAgent(agentName, context, projectId) {
  const response = await apiClient.post(`/agent/${agentName}/query`, {
    context,
    ...(projectId ? { project_id: projectId } : {}),
  }, { timeout: 120000 })
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Project/Campaign/Phase CRUD
// ---------------------------------------------------------------------------

export async function v9ListProjects() {
  const response = await apiClient.get('/v9/projects')
  return response.data
}

export async function v9CreateProject(data) {
  const response = await apiClient.post('/v9/projects', data)
  return response.data
}

export async function v9GetProject(projectId) {
  const response = await apiClient.get(`/v9/projects/${projectId}`)
  return response.data
}

export async function v9UpdateProject(projectId, data) {
  const response = await apiClient.put(`/v9/projects/${projectId}`, data)
  return response.data
}

export async function v9DeleteProject(projectId) {
  await apiClient.delete(`/v9/projects/${projectId}`)
}

export async function v9ListCampaigns(projectId) {
  const response = await apiClient.get(`/v9/projects/${projectId}/campaigns`)
  return response.data
}

export async function v9CreateCampaign(projectId, data) {
  const response = await apiClient.post(`/v9/projects/${projectId}/campaigns`, data)
  return response.data
}

export async function v9GetCampaign(campaignId) {
  const response = await apiClient.get(`/v9/campaigns/${campaignId}`)
  return response.data
}

export async function v9UpdateCampaign(campaignId, data) {
  const response = await apiClient.put(`/v9/campaigns/${campaignId}`, data)
  return response.data
}

export async function v9ListPhases(campaignId) {
  const response = await apiClient.get(`/v9/campaigns/${campaignId}/phases`)
  return response.data
}

export async function v9CreatePhase(campaignId, data) {
  const response = await apiClient.post(`/v9/campaigns/${campaignId}/phases`, data)
  return response.data
}

export async function v9GetPhase(phaseId) {
  const response = await apiClient.get(`/v9/phases/${phaseId}`)
  return response.data
}

export async function v9UpdatePhase(phaseId, data) {
  const response = await apiClient.put(`/v9/phases/${phaseId}`, data)
  return response.data
}

export async function v9FreezePhase(phaseId) {
  const response = await apiClient.post(`/v9/phases/${phaseId}/freeze`)
  return response.data
}

export async function v9UnfreezePhase(phaseId) {
  const response = await apiClient.post(`/v9/phases/${phaseId}/unfreeze`)
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Runs
// ---------------------------------------------------------------------------

export async function v9CreateRun(phaseId, data) {
  const response = await apiClient.post(`/v9/phases/${phaseId}/runs`, data)
  return response.data
}

export async function v9ListRuns(phaseId) {
  const response = await apiClient.get(`/v9/phases/${phaseId}/runs`)
  return response.data
}

export async function v9GetRun(runId) {
  const response = await apiClient.get(`/v9/runs/${runId}`)
  return response.data
}

export async function v9CancelRun(runId) {
  const response = await apiClient.post(`/v9/runs/${runId}/cancel`)
  return response.data
}

export async function v9ArchiveRun(runId) {
  const response = await apiClient.post(`/v9/runs/${runId}/archive`)
  return response.data
}

export async function v9ImportFile(phaseId, file) {
  const formData = new FormData()
  formData.append('file', file)
  const response = await apiClient.post(
    `/v9/phases/${phaseId}/runs/import-file`,
    formData,
    { headers: { 'Content-Type': 'multipart/form-data' }, timeout: 120000 }
  )
  return response.data
}

export async function v9ImportDatabase(phaseId, config) {
  const response = await apiClient.post(
    `/v9/phases/${phaseId}/runs/import-database`,
    config,
    { timeout: 120000 }
  )
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Molecules
// ---------------------------------------------------------------------------

export async function v9ListMolecules(phaseId, params = {}) {
  const response = await apiClient.get(`/v9/phases/${phaseId}/molecules`, { params })
  return response.data
}

export async function v9MoleculeStats(phaseId) {
  const response = await apiClient.get(`/v9/phases/${phaseId}/molecules/stats`)
  return response.data
}

export async function v9GetMolecule(moleculeId) {
  const response = await apiClient.get(`/v9/molecules/${moleculeId}`)
  return response.data
}

export async function v9BookmarkMolecule(moleculeId, bookmarked) {
  const response = await apiClient.put(
    `/v9/molecules/${moleculeId}/bookmark`,
    null,
    { params: { bookmarked } }
  )
  return response.data
}

export async function v9BookmarkBatch(phaseId, moleculeIds, bookmarked) {
  const response = await apiClient.post(`/v9/phases/${phaseId}/molecules/bookmark-batch`, {
    molecule_ids: moleculeIds,
    bookmarked,
  })
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Utilities
// ---------------------------------------------------------------------------

export async function v9SmilesToMolblock(smiles) {
  const response = await apiClient.post('/v9/molecule/smiles-to-3d', { smiles })
  return response.data
}

export default apiClient
