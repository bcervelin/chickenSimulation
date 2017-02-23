##Copyright 2017 Bruno Henrique Cervelin
##This file is part of chickenSimulation
##
##chickenSimulation is free software: you can redistribute it and/or modify
##it under the terms of the GNU General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##This program is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU General Public License for more details.
##
##You should have received a copy of the GNU General Public License
##along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
.SUFFIXES:
.SUFFIXES:  .f90 .o

COMPILER ?= gfortran

simulacao: ff.o aleatory.o update_parameters.o forces.o simu_fortran.o
	$(COMPILER) $^ \
	-o simulacao
#implicit build rule
.f90.o:
	$(COMPILER) -c $^
#cleaning
clean:
	-rm *.o *.mod simulacao *.m *.out
cleaner: clean
cleanest: cleaner
## EXPLICIT DEPENDENCIES
