from pathlib import Path
import html
import sys
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.platypus import Paragraph, Table, TableStyle
from reportlab.lib.enums import TA_LEFT

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'pass to VD-docs/Mass estimator equation'
FONT_DIR = Path('C:/Windows/Fonts')
for name, file in [('Body', 'calibri.ttf'), ('BodyBold', 'calibrib.ttf'),
                   ('BodyItalic', 'calibrii.ttf'), ('Math', 'times.ttf'),
                   ('MathItalic', 'timesi.ttf'), ('MathBold', 'timesbd.ttf'),
                   ('Mono', 'consola.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONT_DIR/file)))
pdfmetrics.registerFontFamily('Body', normal='Body', bold='BodyBold', italic='BodyItalic', boldItalic='BodyBold')
NAVY = '#162B3D'
TEAL = '#087F80'
GREY = '#50616F'
PALE = '#F1F6F7'
LIGHT = '#D8E3E6'
SCENE_W, SCENE_H = 1400, 700
scene = []


def text(x, y, value, size=40, font='Math', color=NAVY):
    scene.append(('text', x, y, value, size, font, color))


def line(x1, y1, x2, y2, width=2, color=NAVY):
    scene.append(('line', x1, y1, x2, y2, width, color))


def hat(x, y, width=29):
    line(x, y, x+width/2, y-8, 2)
    line(x+width/2, y-8, x+width, y, 2)


def integral_row(x, baseline, variable):
    text(x, baseline+7, '∫', 65)
    text(x+22, baseline+19, 'C', 26, 'MathItalic')
    text(x+65, baseline, variable, 50, 'MathItalic')
    text(x+119, baseline, 'e', 50, 'MathItalic')
    text(x+145, baseline-27, '−Q/2', 29, 'MathItalic')
    text(x+240, baseline, 'dν', 46, 'MathItalic')


def bracket(x, y, height=76, right=False):
    edge = x-10 if right else x+10
    line(x, y, edge, y, 1.5)
    line(x, y, x, y+height, 1.5)
    line(x, y+height, edge, y+height, 1.5)


# One scene produces both the editable SVG and vector drawing in the PDF.
text(46, 52, 'THE COMPACT EQUATION  /  GAUSSIAN ERRORS + FLAT POSITIVE MAGNITUDES', 23, 'BodyBold', TEAL)
text(136, 206, 'm', 66, 'MathItalic'); hat(140, 157, 35)
text(210, 206, '=', 53)
text(296, 157, 'E', 51); text(330, 168, 'P', 26, 'MathItalic'); text(352, 157, '[f]', 51)
line(280, 180, 419, 180, 2.3)
text(291, 235, 'E', 51); text(325, 246, 'P', 26, 'MathItalic'); text(347, 235, '[α]', 51)
text(466, 206, '=', 53)
integral_row(574, 149, 'f')
line(555, 180, 936, 180, 2.3)
integral_row(574, 238, 'α')
text(995, 162, 'Average the', 28, 'Body')
text(995, 199, 'magnitudes,', 28, 'Body')
text(995, 236, 'then divide.', 28, 'BodyBold', TEAL)
line(46, 286, 1354, 286, 1.2, LIGHT)

text(64, 339, 'C = (0, ∞)', 37, 'MathItalic'); text(236, 320, '2', 23)
text(265, 339, '× S', 37); text(333, 320, 'd−1', 23, 'MathItalic')
text(492, 339, 'dν = df dα dΩ(u)', 37, 'MathItalic')
text(1000, 339, 'd = 1, 2 or 3', 34, 'MathItalic')
text(64, 406, 'dP = Z', 37, 'MathItalic'); text(183, 386, '−1', 23)
text(224, 406, 'e', 37, 'MathItalic'); text(245, 386, '−Q/2', 23, 'MathItalic')
text(321, 406, 'dν', 37, 'MathItalic')
text(493, 406, 'Z =', 37, 'MathItalic'); text(565, 411, '∫', 46)
text(585, 419, 'C', 22, 'MathItalic'); text(622, 406, 'e', 37, 'MathItalic')
text(644, 386, '−Q/2', 23, 'MathItalic'); text(725, 406, 'dν', 37, 'MathItalic')
text(64, 478, 'Q = (y − x)', 40, 'MathItalic'); text(260, 455, 'T', 24)
text(288, 478, 'Σ', 40); text(317, 455, '−1', 24)
text(354, 478, '(y − x)', 40, 'MathItalic')
text(644, 478, 'Σ: symmetric positive-definite covariance', 26, 'Body')

text(67, 573, 'y =', 39, 'MathItalic')
bracket(145, 522); bracket(226, 522, right=True)
text(169, 550, 'F', 34, 'MathBold'); text(169, 526, '~', 26)
text(171, 587, 'a', 34, 'MathBold'); text(169, 563, '~', 26)
text(351, 573, 'x =', 39, 'MathItalic')
bracket(432, 522); bracket(537, 522, right=True)
text(455, 549, 'fu', 34, 'MathItalic')
text(451, 586, 'αu', 34, 'MathItalic')
text(661, 548, 'f > 0, α > 0, ||u|| = 1', 34, 'MathItalic')
text(661, 589, 'Uniform direction; flat in f and α.', 29, 'Body')
line(46, 624, 1354, 624, 1.2, LIGHT)
text(64, 670, 'P is the normalized law over compatible latent pairs.  F and a above are measured vectors.', 27, 'Body', GREY)


def export_svg():
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{SCENE_W}" height="{SCENE_H}" viewBox="0 0 {SCENE_W} {SCENE_H}" role="img" aria-label="Mass equals the ratio of mean compatible force and acceleration magnitudes; Gaussian flat-magnitude law fully specified below.">',
             '<rect width="100%" height="100%" fill="white"/>']
    for item in scene:
        if item[0] == 'text':
            _, x, y, value, size, font, color = item
            family = 'Calibri, Arial, sans-serif' if font.startswith('Body') else 'Times New Roman, serif'
            italic = 'italic' if 'Italic' in font else 'normal'
            weight = '700' if 'Bold' in font else '400'
            parts.append(f'<text x="{x}" y="{y}" font-family="{family}" font-size="{size}" font-style="{italic}" font-weight="{weight}" fill="{color}">{html.escape(value)}</text>')
        else:
            _, x1, y1, x2, y2, width, color = item
            parts.append(f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{color}" stroke-width="{width}"/>')
    parts.append('</svg>')
    (OUT/'equation.svg').write_text('\n'.join(parts), encoding='utf-8')


def draw_scene(c, x, top, width):
    scale = width/SCENE_W
    c.saveState(); c.translate(x, top-SCENE_H*scale); c.scale(scale, scale)
    for item in scene:
        if item[0] == 'text':
            _, tx, ty, value, size, font, color = item
            c.setFillColor(colors.HexColor(color)); c.setFont(font, size)
            c.drawString(tx, SCENE_H-ty, value)
        else:
            _, x1, y1, x2, y2, lw, color = item
            c.setStrokeColor(colors.HexColor(color)); c.setLineWidth(lw)
            c.line(x1, SCENE_H-y1, x2, SCENE_H-y2)
    c.restoreState()


PAGE_W, PAGE_H = 595.276, 841.89
MARGIN, WIDTH = 48, PAGE_W-96
style = ParagraphStyle('body', fontName='Body', fontSize=11, leading=15,
                       textColor=colors.HexColor(NAVY), spaceAfter=7)
small_style = ParagraphStyle('small', parent=style, fontSize=9.3, leading=12.3)


def para(c, text_value, x, top, width=WIDTH, small=False, size=None):
    st = small_style if small else style
    if size:
        st = ParagraphStyle('custom', parent=st, fontSize=size, leading=size*1.35)
    p = Paragraph(text_value, st)
    w, h = p.wrap(width, PAGE_H)
    p.drawOn(c, x, top-h)
    return top-h-8


def heading(c, title, top):
    c.setFillColor(colors.HexColor(TEAL)); c.setFont('BodyBold', 13)
    c.drawString(MARGIN, top, title)
    return top-15


def base(c, page, title, subtitle):
    c.setFillColor(colors.HexColor(TEAL)); c.rect(MARGIN, PAGE_H-49, 27, 4, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(GREY)); c.setFont('BodyBold', 9)
    c.drawString(MARGIN+39, PAGE_H-49, 'VD NEWTON  /  MASS ESTIMATOR')
    c.setFillColor(colors.HexColor(NAVY)); c.setFont('BodyBold', 27)
    c.drawString(MARGIN, PAGE_H-96, title)
    para(c, subtitle, MARGIN, PAGE_H-113, size=11.3)
    c.setStrokeColor(colors.HexColor(LIGHT)); c.setLineWidth(.6)
    c.line(MARGIN, 45, PAGE_W-MARGIN, 45)
    c.setFillColor(colors.HexColor(GREY)); c.setFont('Body', 8.5)
    c.drawString(MARGIN, 29, 'Formula specification  |  21 September 2026')
    c.drawRightString(PAGE_W-MARGIN, 29, f'{page} / 3')


def table(c, rows, top, widths):
    cellstyle = ParagraphStyle('cell', parent=small_style, fontSize=10, leading=13)
    data = [[Paragraph(cell, cellstyle) for cell in row] for row in rows]
    t = Table(data, colWidths=widths, hAlign='LEFT')
    t.setStyle(TableStyle([
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor(PALE)),
        ('LINEBELOW', (0,0), (-1,0), .7, colors.HexColor(LIGHT)),
        ('LINEBELOW', (0,1), (-1,-1), .4, colors.HexColor(LIGHT)),
        ('LEFTPADDING', (0,0), (-1,-1), 8), ('RIGHTPADDING', (0,0), (-1,-1), 8),
        ('TOPPADDING', (0,0), (-1,-1), 7), ('BOTTOMPADDING', (0,0), (-1,-1), 7)
    ]))
    _, h = t.wrap(WIDTH, PAGE_H)
    t.drawOn(c, MARGIN, top-h)
    return top-h-18


def build_pdf():
    c = canvas.Canvas(str(OUT/'mass_estimator.pdf'), pagesize=(PAGE_W, PAGE_H))
    c.setTitle('The mass estimator equation')
    c.setAuthor('VD Newton')
    c.setSubject('Compact specification and Python implementation of the compatible-pair ratio-of-means mass estimator')
    base(c, 1, 'The mass estimator equation',
         'A positive inertial-mass estimate from measured force and acceleration, with uncertainty in both vectors.')
    draw_scene(c, MARGIN-4, 660, WIDTH+8)
    y = 386
    y = heading(c, 'Types and form', y)
    y = para(c, 'The inputs are two <b>vectors</b> and their joint covariance <b>matrix</b>. The output is a positive <b>scalar mass</b>, accompanied by a conditional interval and a logarithmic spread. The code supports one, two or three spatial dimensions.', MARGIN, y)
    y = para(c, 'The integration variables <font face="MathItalic">f</font> and <font face="Math">α</font> are possible true magnitudes; <font face="MathItalic">u</font> is their common unit direction. The formula averages these magnitudes under one joint probability law, then divides. Its mathematical form is a ratio of weighted integrals.', MARGIN, y)
    y = heading(c, 'Basic physical assumption', y-10)
    y = para(c, '<b>Newton II holds for the true quantities:</b> <font face="MathItalic">F*</font> = <font face="MathItalic">m a*</font>, with <font face="MathItalic">m</font> &gt; 0. Both measurements refer to the same object, time or interval, and inertial frame. Neither input may already assume the mass being sought.', MARGIN, y)
    y = heading(c, 'Statistical assumptions that complete the equation', y-10)
    y = para(c, 'Measurement errors are jointly Gaussian with zero mean and the supplied positive-definite covariance. The common direction is uniform, and the reference measure is flat in the two positive physical magnitudes. These choices are explicit parts of the estimator; Newton II alone does not determine them.', MARGIN, y)
    assert y > 50, y
    c.showPage()

    base(c, 2, 'Inputs and uncertainty',
         'Preserve the vector components and specify their errors in the same coordinate system.')
    rows = [
        ['<b>Input</b>', '<b>Meaning and units</b>'],
        ['<b>force</b><br/>shape (d,)', 'Measured <b>net force</b> components. SI unit: N. Zero and negative components are valid.'],
        ['<b>acceleration</b><br/>shape (d,)', 'Measured acceleration components in the same axes. SI unit: m/s².'],
        ['<b>covariance</b><br/>shape (2d, 2d)', 'Error covariance in <b>[F1, ..., Fd, a1, ..., ad]</b> order. Diagonal entries are variances. Supply cross-correlations when present.'],
        ['<b>probability</b><br/>default 0.95', 'Probability content of the central equal-tail interval. Must lie strictly between zero and one.']
    ]
    y = table(c, rows, 672, [128, WIDTH-128])
    y = para(c, 'For independent isotropic 3D errors, use diagonal entries [σ<sub>F</sub>², σ<sub>F</sub>², σ<sub>F</sub>², σ<sub>a</sub>², σ<sub>a</sub>², σ<sub>a</sub>²]. Each σ is a <b>per-coordinate standard uncertainty</b>. The covariance block units are N², (m/s²)² and N·m/s².', MARGIN, y)
    y = heading(c, 'Uncertainty comes from the same joint law', y-7)
    y = para(c, 'Let <font face="MathItalic">M = f/α</font> under <font face="MathItalic">P</font>, and let <font face="MathItalic">q<sub>p</sub></font> be its probability-<font face="MathItalic">p</font> quantile. The returned uncertainty is:', MARGIN, y)
    c.setFillColor(colors.HexColor(PALE)); c.roundRect(MARGIN, y-65, WIDTH, 59, 6, fill=1, stroke=0)
    para(c, '<b>95% conditional interval:</b> [q<sub>0.025</sub>, q<sub>0.975</sub>]<br/><b>Logarithmic spread:</b> SD<sub>P</sub>[log(M/m<sub>0</sub>)]', MARGIN+14, y-17, WIDTH-28, size=12)
    y -= 80
    y = para(c, 'The reference mass m<sub>0</sub> makes the logarithm dimensionless; its choice does not change the spread. A 95% interval contains 95% of this declared conditional law. It does <b>not automatically promise 95% coverage</b> across repeated experiments.', MARGIN, y)
    y = para(c, 'The ordinary mean and second moment of M diverge in this model, so an ordinary mass standard deviation is undefined. Quantiles and log spread remain finite. Report the interval as endpoints, not as a symmetric mass ± SD.', MARGIN, y)
    y = heading(c, 'Use and interpretation', y-7)
    y = para(c, '<b>1.</b> Establish the object, frame, interval, units and covariance.<br/><b>2.</b> Supply the measured vector components directly.<br/><b>3.</b> Report the point, interval, probability and declared model together.<br/><b>4.</b> Keep the original measurements and covariance with their provenance.', MARGIN, y)
    y = para(c, 'The point need not equal the median or interval midpoint. Inverse mass uses 1/m̂ and the reversed reciprocal interval [1/upper, 1/lower].', MARGIN, y, small=True)
    assert y > 50, y
    c.showPage()

    base(c, 3, 'Calculate and report it',
         'A standalone Python function is included in this folder. It requires Python 3.10+ and NumPy.')
    y = 670
    y = para(c, '<font face="BodyBold">Install:</font> <font face="Mono">python -m pip install -r requirements.txt</font><br/><font face="BodyBold">Run the supplied example:</font> <font face="Mono">python example.py</font>', MARGIN, y, small=True)
    code = [
        'import numpy as np',
        'from mass_estimator import estimate_mass',
        '',
        'force = [10.0, 0.0, 0.0]        # N',
        'acceleration = [2.0, 0.1, 0.0] # m/s^2',
        'cov = np.diag([0.5**2]*3 + [0.1**2]*3)',
        'r = estimate_mass(force, acceleration, cov)',
        'print(r.mass, r.interval)',
        'print(r.log_standard_deviation)',
    ]
    height = 145
    c.setFillColor(colors.HexColor(PALE)); c.roundRect(MARGIN, y-height, WIDTH, height, 6, fill=1, stroke=0)
    c.setFillColor(colors.HexColor(NAVY)); c.setFont('Mono', 9.1)
    for i, row in enumerate(code): c.drawString(MARGIN+14, y-21-i*13.5, row)
    y -= height+18
    y = table(c, [
        ['<b>Output</b>', '<b>Example value</b>'],
        ['Mass point', '<b>4.993754 kg</b>'],
        ['95% conditional interval', '<b>[4.344146, 5.740363] kg</b>'],
        ['Log standard deviation', '<b>0.071035</b> (dimensionless)']
    ], y, [210, WIDTH-210])
    y = heading(c, 'Scope and numerical checks', y)
    y = para(c, 'The function accepts full covariance, including anisotropy and correlated channels. An optional externally known direction replaces the direction average. The supplied covariance is treated as fixed; uncertainty in its calibration requires the separate extended model.', MARGIN, y, small=True)
    y = para(c, 'Measured zeros and nonaligned vectors are valid inputs. Exact quantities, hard bounds and changing-excitation joint trials require separate treatment. At true zero excitation, measurements cannot identify mass; a finite model-based output does not establish identification.', MARGIN, y, small=True)
    y = para(c, 'Every returned result passes a grid-refinement comparison. The example’s final scaled change is 6.61 × 10<super>−5</super>; this is a numerical diagnostic, not mass uncertainty or a certified error bound. Unresolved integration raises an error. Nine verification checks passed, including analytic identities and independent integrations.', MARGIN, y, small=True)
    y = heading(c, 'Which version this specifies', y-5)
    y = para(c, 'The current <b>ratio of magnitude means</b> is retained in the local 21 September synthesis. The earlier best-point and median proposals were checked. The 7 September steering assumes Newton II and removes empirical rejection from this conditional use.', MARGIN, y, small=True)
    sha = '8782a293297330f7a354d87ebe62792018fa8866'
    prefix = f'https://github.com/canvalidk/VD-docs/blob/{sha}/Newton-analysis/_dscn_mass_estimator/'
    sources = (f'Managed sources: <a href="{prefix}inertial_mass_estimator_full_working_spec_d1.md" color="{TEAL}">working specification §§7-12, 22-23</a>; '
               f'<a href="{prefix}mass_estimator_law_assumed_correct_steering_2026-09-07.md" color="{TEAL}">7 September steering</a>; '
               f'<a href="{prefix}what_forces_the_ratio_of_means_readout.md" color="{TEAL}">readout audit</a>. '
               f'VD-docs commit: {sha}. Full source paths, local source hashes, version distinctions and the API contract are in README.md.')
    y = para(c, sources, MARGIN, y, small=True)
    y = para(c, '<b>Status:</b> user-requested presentation, staged in pass to VD-docs; awaiting transfer.', MARGIN, y, small=True)
    assert y > 50, y
    c.save()


export_svg()
build_pdf()
print(OUT/'mass_estimator.pdf')
