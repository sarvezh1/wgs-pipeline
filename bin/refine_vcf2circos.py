#!/usr/bin/env python3

import argparse
import json
import math
from pathlib import Path

from plotly.offline import plot


def plotly_arguments(html):
    start = html.rfind('Plotly.newPlot(') + len('Plotly.newPlot(')
    end = html.find(')                };', start)
    if start < len('Plotly.newPlot(') or end < 0:
        raise ValueError('Plotly figure payload not found')

    payload = html[start:end]
    parts = []
    token_start = 0
    depth = 0
    quoted = False
    escaped = False
    for index, character in enumerate(payload):
        if quoted:
            if escaped:
                escaped = False
            elif character == '\\':
                escaped = True
            elif character == '"':
                quoted = False
        else:
            if character == '"':
                quoted = True
            elif character in '[{(':
                depth += 1
            elif character in ']})':
                depth -= 1
            elif character == ',' and depth == 0:
                parts.append(payload[token_start:index].strip())
                token_start = index + 1
    parts.append(payload[token_start:].strip())
    return json.loads(parts[1]), json.loads(parts[2])


parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--sample', required=True)
parser.add_argument('--caller', required=True)
parser.add_argument('--analysis-type', required=True)
parser.add_argument('--selection', required=True)
parser.add_argument('--chromosome', default='chr20')
parser.add_argument('--start', type=int, required=True)
parser.add_argument('--end', type=int, required=True)
parser.add_argument('--chromosome-length', type=int, required=True)
parser.add_argument('--records', type=int, required=True)
parser.add_argument('--snvs', type=int, required=True)
parser.add_argument('--indels', type=int, required=True)
args = parser.parse_args()

data, layout = plotly_arguments(Path(args.input).read_text())

trace_counts = {}
original_angles = {}
for trace in data:
    name = trace.get('name')
    if name in {'SNV', 'INDEL'}:
        trace_counts[name] = len(trace.get('text', [])) // 2
        original_angles[name] = [math.atan2(y, x) for x, y in zip(trace['x'], trace['y'])]

        target_radius = 0.84 if name == 'SNV' else 0.91
        for index, (x_value, y_value) in enumerate(zip(trace['x'], trace['y'])):
            scale = target_radius / math.hypot(x_value, y_value)
            trace['x'][index] = x_value * scale
            trace['y'][index] = y_value * scale

        trace['marker']['size'] = 12 if name == 'SNV' else 15
        trace['marker']['symbol'] = 'circle' if name == 'SNV' else 'diamond'
        trace['marker']['color'] = '#356FA8' if name == 'SNV' else '#D97706'
        trace['marker']['line'] = {'color': 'white', 'width': 1}
        trace['legendrank'] = 1 if name == 'SNV' else 2
        trace['showlegend'] = True
    elif name in {'INTERMEDIATE', 'CYTOBAND'}:
        trace['showlegend'] = False

expected_counts = {'SNV': args.snvs, 'INDEL': args.indels}
if trace_counts != expected_counts or sum(trace_counts.values()) != args.records:
    raise ValueError(f'Unexpected variant trace counts: {trace_counts}; expected {expected_counts}')

for trace in data:
    name = trace.get('name')
    if name in original_angles:
        updated_angles = [math.atan2(y, x) for x, y in zip(trace['x'], trace['y'])]
        if any(abs(before - after) > 1e-12 for before, after in zip(original_angles[name], updated_angles)):
            raise ValueError(f'Angular genomic coordinates changed for {name}')

ring_styles = {
    ('lightgray', 'lightgray', 1),
    ('white', 'white', 1),
    ('grey', 'grey', 1),
}
layout['shapes'] = [
    shape for shape in layout.get('shapes', [])
    if (shape.get('fillcolor'), shape.get('line', {}).get('color'), shape.get('line', {}).get('width'))
    not in ring_styles
]

spacing = 0.05 * args.chromosome_length
total_span = args.chromosome_length + spacing
chromosome_offset = spacing / 2
theta_start = 2 * math.pi * (args.start + chromosome_offset) / total_span
theta_end = 2 * math.pi * (args.end + chromosome_offset) / total_span
theta_values = [theta_start + index * (theta_end - theta_start) / 63 for index in range(64)]
interval_radius = 1.13
interval_x = [math.sin(theta) * interval_radius for theta in theta_values]
interval_y = [math.cos(theta) * interval_radius for theta in theta_values]
interval_label = f'{args.chromosome}:{args.start:,}–{args.end:,}'
interval_kb = (args.end - args.start) / 1000

data.append({
    'type': 'scatter',
    'mode': 'lines',
    'x': interval_x,
    'y': interval_y,
    'name': 'Tested interval',
    'legendrank': 3,
    'showlegend': True,
    'hoverinfo': 'text',
    'text': [f'<b>Technical interval</b><br>{interval_label}<br>GRCh38'] * len(interval_x),
    'line': {'color': '#178F8A', 'width': 14},
})

annotations = []
for annotation in layout.get('annotations', []):
    if annotation.get('text') in {'0 Mb', '30 Mb', '60 Mb'}:
        annotation['font'] = {'family': 'Arial', 'size': 11, 'color': '#4B5563'}
        annotations.append(annotation)

theta_midpoint = (theta_start + theta_end) / 2
highlight_x = math.sin(theta_midpoint) * interval_radius
highlight_y = math.cos(theta_midpoint) * interval_radius
annotations.extend([
    {
        'xref': 'paper', 'yref': 'paper', 'x': 0.5, 'y': 0.985,
        'text': '<b>GRCh38 · full chromosome 20 ideogram</b>',
        'showarrow': False,
        'font': {'family': 'Arial', 'size': 15, 'color': '#374151'},
    },
    {
        'xref': 'x', 'yref': 'y', 'x': 0, 'y': 0.08,
        'text': f'<b>{args.chromosome}</b>', 'showarrow': False,
        'font': {'family': 'Arial', 'size': 34, 'color': '#1F2937'},
    },
    {
        'xref': 'x', 'yref': 'y', 'x': 0, 'y': -0.04,
        'text': f'<b>{interval_kb:.1f} kb tested interval</b><br>{args.snvs} SNVs · {args.indels} indels',
        'showarrow': False,
        'font': {'family': 'Arial', 'size': 17, 'color': '#4B5563'},
    },
    {
        'xref': 'paper', 'yref': 'paper', 'x': 0.5, 'y': 0.015,
        'text': f'Highlighted arc marks the {interval_kb:.1f} kb technical interval · {interval_label}',
        'showarrow': False,
        'font': {'family': 'Arial', 'size': 14, 'color': '#374151'},
        'bgcolor': 'rgba(255,255,255,0.92)',
        'bordercolor': '#D1D5DB', 'borderwidth': 1, 'borderpad': 6,
    },
    {
        'xref': 'x', 'yref': 'y', 'x': highlight_x, 'y': highlight_y,
        'text': '<b>tested interval</b>', 'showarrow': True,
        'arrowhead': 2, 'arrowwidth': 2, 'arrowcolor': '#178F8A',
        'ax': 72, 'ay': 48,
        'font': {'family': 'Arial', 'size': 13, 'color': '#176B68'},
        'bgcolor': 'rgba(255,255,255,0.9)', 'borderpad': 3,
    },
])

selection_label = 'PASS variants' if args.selection == 'PASS' else 'variants'
interval_mb = f'{args.start / 1_000_000:.2f}–{args.end / 1_000_000:.2f} Mb'
layout.update({
    'width': 1200,
    'height': 1000,
    'paper_bgcolor': 'white',
    'plot_bgcolor': 'white',
    'font': {'family': 'Arial', 'color': '#1F2937'},
    'margin': {'l': 45, 'r': 45, 't': 120, 'b': 90},
    'annotations': annotations,
    'title': {
        'text': (
            f'<b>{args.sample} {args.analysis_type} variants | {args.caller.capitalize()}</b>'
            f'<br><span style="font-size:17px;color:#4B5563">'
            f'GRCh38 {args.chromosome}:{interval_mb} | {args.records} {selection_label}</span>'
        ),
        'x': 0.5, 'xanchor': 'center', 'y': 0.975, 'yanchor': 'top',
        'font': {'family': 'Arial', 'size': 25, 'color': '#111827'},
    },
    'legend': {
        'title': {'text': '<b>Variant class</b>', 'font': {'family': 'Arial', 'size': 15}},
        'font': {'family': 'Arial', 'size': 15, 'color': '#1F2937'},
        'orientation': 'v', 'x': 0.985, 'xanchor': 'right', 'y': 0.94, 'yanchor': 'top',
        'bgcolor': 'rgba(255,255,255,0.92)', 'bordercolor': '#D1D5DB', 'borderwidth': 1,
        'itemsizing': 'constant',
    },
})
layout['xaxis']['range'] = [-1.30, 1.30]
layout['yaxis']['range'] = [-1.30, 1.30]

plot(
    {'data': data, 'layout': layout},
    filename=args.output,
    auto_open=False,
    include_plotlyjs=True,
    config={'responsive': True, 'displaylogo': False, 'displayModeBar': True, 'scrollZoom': True},
)

output_html = Path(args.output).read_text()
for forbidden in ('/tmp/', '/work/', 'file:///'):
    if forbidden in output_html:
        raise ValueError(f'Non-portable path detected: {forbidden}')
