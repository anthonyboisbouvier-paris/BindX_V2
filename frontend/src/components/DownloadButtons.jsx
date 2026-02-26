import React, { useState } from 'react'
import { getReportUrl, getDownloadUrl } from '../api.js'

function DownloadButton({ href, label, icon, variant = 'blue', description }) {
  const [clicking, setClicking] = useState(false)

  const handleClick = () => {
    setClicking(true)
    setTimeout(() => setClicking(false), 1500)
  }

  const baseClasses = 'flex items-center gap-3 px-5 py-3 rounded-xl font-semibold text-sm transition-all duration-200 shadow-sm hover:shadow-md'
  const variantClasses = variant === 'green'
    ? 'bg-dockit-green hover:bg-dockit-green-dark text-white'
    : 'bg-dockit-blue hover:bg-dockit-blue-light text-white'

  return (
    <a
      href={href}
      target="_blank"
      rel="noopener noreferrer"
      download
      onClick={handleClick}
      className={`${baseClasses} ${variantClasses} ${clicking ? 'scale-95' : 'hover:scale-105'}`}
    >
      <div className="flex-shrink-0">
        {icon}
      </div>
      <div className="text-left">
        <div>{clicking ? 'Downloading...' : label}</div>
        {description && (
          <div className="text-xs opacity-75 font-normal">{description}</div>
        )}
      </div>
    </a>
  )
}

const PdfIcon = () => (
  <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
      d="M7 21h10a2 2 0 002-2V9.414a1 1 0 00-.293-.707l-5.414-5.414A1 1 0 0012.586 3H7a2 2 0 00-2 2v14a2 2 0 002 2z" />
    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M9 13h6m-6 4h6m2-10v3a1 1 0 01-1 1h-3" />
    <text x="8" y="19" style={{ fontSize: '6px', fontWeight: 'bold', fill: 'currentColor' }}>PDF</text>
  </svg>
)

const ZipIcon = () => (
  <svg className="w-6 h-6" fill="none" stroke="currentColor" viewBox="0 0 24 24">
    <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
      d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-4l-4 4m0 0l-4-4m4 4V4" />
  </svg>
)

export default function DownloadButtons({ jobId, results }) {
  const reportUrl = getReportUrl(jobId)
  const downloadUrl = getDownloadUrl(jobId)

  const totalResults = results?.results?.length || 0
  const jobInfo = results?.uniprot_id || jobId

  return (
    <div className="card p-5">
      {/* Header */}
      <div className="flex items-center gap-2 mb-4">
        <svg className="w-5 h-5 text-dockit-blue" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2}
            d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z" />
        </svg>
        <h3 className="font-semibold text-dockit-blue">Download Results</h3>
      </div>

      {/* Summary info */}
      <div className="flex items-center gap-4 mb-4 p-3 bg-green-50 rounded-lg border border-green-100">
        <div className="w-8 h-8 bg-dockit-green rounded-full flex items-center justify-center flex-shrink-0">
          <svg className="w-4 h-4 text-white" fill="currentColor" viewBox="0 0 20 20">
            <path fillRule="evenodd" d="M16.707 5.293a1 1 0 010 1.414l-8 8a1 1 0 01-1.414 0l-4-4a1 1 0 011.414-1.414L8 12.586l7.293-7.293a1 1 0 011.414 0z" clipRule="evenodd" />
          </svg>
        </div>
        <div>
          <p className="text-sm font-semibold text-green-800">
            Screening completed successfully
          </p>
          <p className="text-xs text-green-600">
            {totalResults} ligands tested for {jobInfo}
          </p>
        </div>
      </div>

      {/* Buttons */}
      <div className="flex flex-col sm:flex-row gap-3">
        <DownloadButton
          href={reportUrl}
          label="PDF Report"
          description="Full summary with figures"
          icon={<PdfIcon />}
          variant="blue"
        />
        <DownloadButton
          href={downloadUrl}
          label="All Results (ZIP)"
          description="PDB, PDBQT, CSV, report"
          icon={<ZipIcon />}
          variant="green"
        />
      </div>

      {/* Contents description */}
      <div className="mt-4 grid grid-cols-2 gap-2">
        <div className="p-2.5 bg-gray-50 rounded-lg">
          <p className="text-xs font-semibold text-gray-600 mb-1">The PDF report contains:</p>
          <ul className="text-xs text-gray-400 space-y-0.5">
            <li>• Screening summary</li>
            <li>• Top 10 ligands with scores</li>
            <li>• 2D molecular structures</li>
            <li>• Physicochemical properties</li>
          </ul>
        </div>
        <div className="p-2.5 bg-gray-50 rounded-lg">
          <p className="text-xs font-semibold text-gray-600 mb-1">The ZIP archive contains:</p>
          <ul className="text-xs text-gray-400 space-y-0.5">
            <li>• Receptor PDB structure</li>
            <li>• PDBQT poses for each ligand</li>
            <li>• Complete CSV data</li>
            <li>• PDF report</li>
          </ul>
        </div>
      </div>
    </div>
  )
}
