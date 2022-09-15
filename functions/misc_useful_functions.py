#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:50:42 2022

@author: qinqinyu
"""

import numpy as np
from epiweeks import Week, Year
from datetime import timedelta

def epiweeks2dates(data_epiweeks):
    # Converts array of epiweeks to array of dates
    dates = []
    for w in data_epiweeks:
        w = int(np.floor(w))
        if w<=53:
            week = Week(2020, w)
        elif w<=105:
            week = Week(2021, w-53)
        else:
            week = Week(2022, w-105)
        dates.append(week.startdate()+timedelta(days=3))
        
    dates = np.array(dates)
    return dates