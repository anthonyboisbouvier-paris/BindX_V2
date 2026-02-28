import axios from 'axios'

const BASE_URL = import.meta.env.VITE_API_URL || '/api'

const apiClient = axios.create({
  baseURL: BASE_URL,
  timeout: 30000,
  headers: {
    'Content-Type': 'application/json',
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

/**
 * Create a new docking job
 * @param {Object} params - Job parameters
 * @param {string} params.uniprot_id - UniProt ID (e.g. P00533)
 * @param {string} [params.sequence] - Raw amino acid sequence (alternative to uniprot_id)
 * @param {string[]} [params.custom_smiles] - Custom SMILES strings
 * @param {boolean} [params.use_chembl] - Use ChEMBL ligands
 * @param {boolean} [params.use_zinc] - Use ZINC ligands
 * @param {number} [params.max_ligands] - Max number of ligands
 * @param {string} [params.mode] - Analysis mode: 'basic' | 'advanced'
 * @param {string} [params.docking_engine] - Engine: 'vina' | 'diffdock'
 * @param {string} [params.notification_email] - Optional email for job completion notification
 * @returns {Promise<{job_id: string}>}
 */
export async function createJob(params) {
  // Pass all V3 fields through as-is; the backend ignores unknown fields gracefully
  const response = await apiClient.post('/jobs', params)
  return response.data
}

/**
 * Get job status and progress
 * @param {string} jobId
 * @returns {Promise<{job_id: string, status: string, progress: number, current_step: string, error?: string}>}
 */
export async function getJobStatus(jobId) {
  const response = await apiClient.get(`/jobs/${jobId}`)
  return response.data
}

/**
 * Get job results
 * @param {string} jobId
 * @returns {Promise<Object>} Full results object
 */
export async function getJobResults(jobId) {
  const response = await apiClient.get(`/jobs/${jobId}/results`)
  return response.data
}

/**
 * Get URL for downloading the PDF report
 * @param {string} jobId
 * @returns {string}
 */
export function getReportUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/report`
}

/**
 * Get URL for downloading the ZIP archive
 * @param {string} jobId
 * @returns {string}
 */
export function getDownloadUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/download`
}

/**
 * Get URL for the protein PDB file
 * @param {string} jobId
 * @returns {string}
 */
export function getProteinUrl(jobId) {
  return `${BASE_URL}/jobs/${jobId}/protein`
}

/**
 * Get URL for a specific docking pose PDBQT file
 * @param {string} jobId
 * @param {number} index - Pose index (0-based)
 * @returns {string}
 */
export function getPoseUrl(jobId, index) {
  return `${BASE_URL}/jobs/${jobId}/pose/${index}`
}

/**
 * Check API health
 * @returns {Promise<{status: string}>}
 */
export async function checkHealth() {
  const response = await apiClient.get('/health')
  return response.data
}

/**
 * Preview target information before submitting a job.
 * Returns structure source, pocket data, and ChEMBL stats.
 * @param {string} uniprotId - UniProt accession ID
 * @returns {Promise<Object>} Preview data
 */
export async function previewTarget(uniprotId) {
  const response = await apiClient.post('/preview-target', { uniprot_id: uniprotId }, { timeout: 120000 })
  return response.data
}

export async function previewSequence(sequence) {
  const response = await apiClient.post('/preview-sequence', { sequence }, { timeout: 120000 })
  return response.data
}

/**
 * Detect pockets on any PDB structure URL using P2Rank.
 * @param {string} downloadUrl - URL to a PDB file
 * @param {string} [label] - Label for logging
 * @param {string} [ligandId] - Co-crystal ligand ID (optional)
 * @returns {Promise<Object>} { pockets: [...] }
 */
export async function detectPockets(downloadUrl, label, ligandId) {
  const response = await apiClient.post('/detect-pockets', {
    download_url: downloadUrl,
    label: label || 'structure',
    ligand_id: ligandId || null,
  }, { timeout: 120000 })
  return response.data
}

/**
 * Get retrosynthesis route for a specific molecule in a job
 * @param {string} jobId
 * @param {number} molIndex - 0-based molecule index in results
 * @returns {Promise<Object>} Synthesis route object
 */
export const getSynthesisRoute = (jobId, molIndex) =>
  apiClient.get(`/jobs/${jobId}/synthesis/${molIndex}`).then(r => r.data)

/**
 * Analyze a molecule's scaffold and R-group positions.
 * @param {string} smiles - SMILES string of the molecule
 * @returns {Promise<Object>} Scaffold analysis with positions, SVG, core indices
 */
export async function analyzeScaffold(smiles) {
  const response = await apiClient.post('/molecule/analyze-scaffold', { smiles })
  return response.data
}

/**
 * Start a lead optimization run for a given molecule in a job.
 * @param {string} jobId
 * @param {Object} params - Optimization parameters
 * @param {Object} params.weights - Objective weights (binding_affinity, toxicity, bioavailability, synthesis)
 * @param {number} params.n_iterations - Number of optimization iterations
 * @param {number} params.variants_per_iteration - Variants per iteration
 * @param {string} [params.starting_smiles] - SMILES of the starting molecule
 * @returns {Promise<{opt_id: string}>}
 */
export async function startOptimization(jobId, params) {
  const response = await apiClient.post(`/jobs/${jobId}/optimize`, params)
  return response.data
}

/**
 * Poll the status of a lead optimization run.
 * @param {string} jobId
 * @param {string} optId - Optimization run ID returned by startOptimization
 * @returns {Promise<Object>} Optimization status object with iterations, best_molecule, objectives
 */
export async function getOptimizationStatus(jobId, optId) {
  const response = await apiClient.get(`/jobs/${jobId}/optimization/${optId}`)
  return response.data
}

/**
 * Retrieve the audit log for a job (traceability of all pipeline steps).
 * @param {string} jobId
 * @returns {Promise<Object>} Audit log entries
 */
export async function getAuditLog(jobId) {
  const response = await apiClient.get(`/jobs/${jobId}/audit_log`)
  return response.data
}

// ---------------------------------------------------------------------------
// V7: Auth API
// ---------------------------------------------------------------------------

export async function authRegister(email, username, password) {
  const response = await apiClient.post('/auth/register', { email, username, password })
  return response.data
}

export async function authLogin(email, password) {
  const response = await apiClient.post('/auth/login', { email, password })
  return response.data
}

export async function authMe() {
  const response = await apiClient.get('/auth/me')
  return response.data
}

// ---------------------------------------------------------------------------
// V7: Projects API
// ---------------------------------------------------------------------------

export async function listProjects() {
  const response = await apiClient.get('/projects')
  return response.data
}

export async function createProject(data) {
  const response = await apiClient.post('/projects', data)
  return response.data
}

export async function getProjectDetail(projectId) {
  const response = await apiClient.get(`/projects/${projectId}`)
  return response.data
}

export async function updateProject(projectId, data) {
  const response = await apiClient.put(`/projects/${projectId}`, data)
  return response.data
}

// ---------------------------------------------------------------------------
// BindX: Target Assessment API
// ---------------------------------------------------------------------------

export async function runTargetAssessment(uniprotId, diseaseContext = null, options = {}) {
  const response = await apiClient.post('/target-assessment', {
    uniprot_id: uniprotId,
    disease_context: diseaseContext,
    ...options,
  }, { timeout: 120000 })
  return response.data
}

export async function getTargetAssessment(assessmentId) {
  const response = await apiClient.get(`/target-assessment/${assessmentId}`)
  return response.data
}

export async function queryAgent(agentName, context, projectId) {
  const response = await apiClient.post(`/agent/${agentName}/query`, {
    context,
    ...(projectId ? { project_id: projectId } : {}),
  }, { timeout: 120000 })
  return response.data
}

export async function triggerRunAnalysis(jobId) {
  const response = await apiClient.post(`/jobs/${jobId}/agent-analysis`, {}, { timeout: 120000 })
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Time Estimation API
// ---------------------------------------------------------------------------

export async function estimatePipelineTime(params) {
  const response = await apiClient.post('/estimate-time', params)
  return response.data
}

export async function listDatabases() {
  const response = await apiClient.get('/databases')
  return response.data
}

// ---------------------------------------------------------------------------
// V9: Project/Campaign/Phase CRUD API
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
// V9: Runs API
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

// ---------------------------------------------------------------------------
// V9: Molecules API
// ---------------------------------------------------------------------------

export async function v9ListMolecules(phaseId, params = {}) {
  const response = await apiClient.get(`/v9/phases/${phaseId}/molecules`, { params })
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

export default apiClient
