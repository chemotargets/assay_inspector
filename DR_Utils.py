#!/usr/bin/env python

__author__ = "Raquel Parrondo-Pizarro"
__date__ = "20250220"
__copyright__ = "Copyright 2024, Chemotargets"
__license__ = ""
__credits__ = ["Data Science & Translational Research Group"]
__maintainer__ = "Raquel Parrondo-Pizarro"
__version__ = "20250220"
__deprecated__ = False

###
def classify_skewness(skewness):
    if skewness < -2.0:
        return 'Severe Left Skewed'
    elif skewness < -1.0:
        return 'Moderate Left Skewed'
    elif skewness > 1.0:
        return 'Moderate Right Skewed'
    elif skewness > 2.0:
        return 'Severe Right Skewed'
    else:
        return 'Non-Skewed'

###
def formatWarningTitle(text):
    base_text = f" {text} Warning "
    num_equals = (75 - len(base_text)) // 2
    equals_str = '=' * num_equals

    # In case the total length is odd, add one more '=' to the end
    if (len(equals_str) * 2 + len(base_text)) < 75:
        return f"\n{equals_str}{base_text}{equals_str}=\n"
    else:
        return f"\n{equals_str}{base_text}{equals_str}\n"