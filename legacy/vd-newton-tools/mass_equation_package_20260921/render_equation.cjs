const path = require('path');
const sharp = require('C:/Users/canva/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/sharp');
const folder = path.resolve(__dirname, '../../pass to VD-docs/Mass estimator equation');
sharp(path.join(folder, 'equation.svg'), {density: 180})
  .png().toFile(path.join(folder, 'equation.png'))
  .then(info => console.log(info.width + 'x' + info.height))
  .catch(error => {console.error(error); process.exitCode = 1;});
