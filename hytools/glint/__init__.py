# -*- coding: utf-8 -*-
"""
HyTools:  Hyperspectral image processing library
Copyright (C) 2021 University of Wisconsin

Authors: Evan Greenberg.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

The :mod:`hytools.correction` module include functions image correction.
"""
from .glint import *
from .gao_2021 import *
from .hedley_2005 import *
from .hochberg_2003 import *
