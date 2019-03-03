# CE2-03-2 Reactor and Controller Design Project

`reactor.m` is the main file. It can be run independently, and will draw constants from `constants.m`(!). `plots.m` does what it says on the tin.

## Notes

- `.slx` files are binary, so `git diff` doesn't work—[use Simulink to compare](https://uk.mathworks.com/help/simulink/ug/merge-simulink-models-from-the-comparison-report.html)
- `functions/` isn't on the path by default—but command `addpath(genpath(pwd))` will add all subfolders

## Useful links

- [`rcdp` on Overleaf](https://www.overleaf.com/project/5c79642b6137e10d5ee8b4ef)
- [`rcdp-equations` on Overleaf](https://www.overleaf.com/project/5c6e80e80219993c39a6dd8a)
- [dynamically access (or create) table / struct data](https://uk.mathworks.com/help/matlab/matlab_prog/generate-field-names-from-variables.html)
