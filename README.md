# CE2-03-2 Reactor and Controller Design Project

### Reactor

- `reactor.m` is the main file. It can be run independently, and will draw constants from `constants.m`(!)
- `plots.m` does what it says on the tin
- `sensitivity.m` tests effect of varying key parameters on the output of `reactor.m`

### Controller

- `PA_Reactor_2019.slx` is the main file
- `SimulinkTest.m` sets up a number of test scenarios

## Notes

- `.slx` files are binary, so `git diff` doesn't work—[use Simulink to compare](https://uk.mathworks.com/help/simulink/ug/merge-simulink-models-from-the-comparison-report.html)
- `functions/` isn't on the path by default—but command `addpath(genpath(pwd))` will add all subfolders

## Useful links

- [`rcdp-report` on Overleaf](https://www.overleaf.com/project/5c79642b6137e10d5ee8b4ef)
- [`rcdp-addendum` on Overleaf](https://www.overleaf.com/project/5c8b7216886bd01953655746)
