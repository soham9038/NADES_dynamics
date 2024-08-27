
## Script to calculate MSD ##

# Software: GROMACS

## Starting of script ##

gmx msd -f nopbc.xtc -s md.tpr -n index.ndx -o msd_glc.xvg
gmx msd -f nopbc.xtc -s md.tpr -n index.ndx -o msd_ure.xvg
gmx msd -f nopbc.xtc -s md.tpr -n index.ndx -o msd_sol.xvg

## End of script ##


