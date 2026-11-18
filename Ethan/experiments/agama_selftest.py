"""Real AGAMA source-adapter smoke test (requires compiled AGAMA)."""
import json
from pathlib import Path
from nbody.galaxies.sources.agama import generate
from nbody.galaxies.agama_status import status
from nbody.io import save

root = Path(__file__).resolve().parents[1]
config = json.loads((root / 'presets/mice.json').read_text())
config['particles'] = {'disk': 300, 'bulge': 100, 'halo': 500}
config['Q'] = 1.5
print('AGAMA environment:', status())
p, report = generate(config, with_report=True)
folder = root / 'results' / 'galaxies'
folder.mkdir(parents=True, exist_ok=True)
save(folder / 'agama_smoke_900.npz', p, {'backend': 'agama', 'config': config,
                                      'generation_report': report})
(folder / 'agama_smoke_900.json').write_text(json.dumps(report, indent=2))
print('Generated:', len(p.m), 'particles, masses:', report['sampled_masses'])
print('DF masses:', report['distribution_function_masses'])
print('Saved to:', folder)
print('IMPORTANT: generation is not stability validation. Run isolation next.')
