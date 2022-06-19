import re

TIER_COLOR_MAPPING = {
    ('143', '25', '27'): "TIER_1",
    ('211', '87', '39'): "TIER_2",
    ('218', '217', '130'): "TIER_3",
    ('253', '251', '250'): "OTHER_MUTATION",
}


def get_tier(style):
    rgb_tuple = re.search(r"rgb\((\d+), (\d+), (\d+)\)", style).groups()
    return TIER_COLOR_MAPPING[rgb_tuple]



