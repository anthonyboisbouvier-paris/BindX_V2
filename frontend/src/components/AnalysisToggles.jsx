import { useRef, useEffect } from 'react';

const TOGGLES = [
  { key: 'enable_admet', label: 'ADMET', desc: 'Absorption, distribution, metabolism, excretion, toxicity' },
  { key: 'enable_synthesis', label: 'Synthesis', desc: 'Retrosynthesis route + cost estimate' },
  { key: 'enable_selectivity', label: 'Selectivity', desc: 'Off-target screening (SEA + structural)' },
  { key: 'enable_herg', label: 'hERG Safety', desc: 'hERG channel inhibition prediction' },
  { key: 'enable_safety', label: 'Full Safety', desc: 'Confidence scoring, PAINS alerts, applicability domain' },
];

const DEFAULT_VALUES = Object.fromEntries(TOGGLES.map(t => [t.key, true]));

export default function AnalysisToggles({ values = DEFAULT_VALUES, onChange }) {
  const masterRef = useRef(null);
  const allChecked = TOGGLES.every(t => values[t.key]);
  const noneChecked = TOGGLES.every(t => !values[t.key]);
  const indeterminate = !allChecked && !noneChecked;

  useEffect(() => {
    if (masterRef.current) {
      masterRef.current.indeterminate = indeterminate;
    }
  }, [indeterminate]);

  const handleMaster = () => {
    const newVal = !allChecked;
    onChange(Object.fromEntries(TOGGLES.map(t => [t.key, newVal])));
  };

  const handleToggle = (key) => {
    onChange({ ...values, [key]: !values[key] });
  };

  return (
    <div className="space-y-2">
      <label className={`flex items-center gap-2 px-3 py-2 rounded-lg border cursor-pointer transition-colors ${
        allChecked ? 'border-[#00e6a0] bg-green-50/40' : 'border-gray-200 bg-white'
      }`}>
        <input
          ref={masterRef}
          type="checkbox"
          checked={allChecked}
          onChange={handleMaster}
          className="w-4 h-4 rounded accent-[#00e6a0]"
        />
        <div>
          <span className="font-medium text-bx-light-text text-sm">All analyses</span>
          <span className="text-xs text-gray-400 ml-2">(recommended)</span>
        </div>
      </label>

      <div className="ml-6 space-y-1">
        {TOGGLES.map(t => (
          <label key={t.key} className={`flex items-start gap-2 px-2 py-1.5 rounded cursor-pointer transition-colors ${
            values[t.key] ? 'bg-green-50/30' : ''
          }`}>
            <input
              type="checkbox"
              checked={!!values[t.key]}
              onChange={() => handleToggle(t.key)}
              className="w-3.5 h-3.5 mt-0.5 rounded accent-[#00e6a0]"
            />
            <div className="leading-tight">
              <span className="text-sm font-medium text-bx-light-text">{t.label}</span>
              <span className="text-xs text-gray-400 ml-1.5">{t.desc}</span>
            </div>
          </label>
        ))}
      </div>
    </div>
  );
}

export { TOGGLES, DEFAULT_VALUES };
