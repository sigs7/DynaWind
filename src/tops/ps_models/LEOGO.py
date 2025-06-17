def load():
    return {
        'base_mva': 100,
        'f': 50,
        'slack_bus': 'Main Bus B',

        'buses': [
            ['name',    'V_n'],
            #['Main Bus A', 11],
            ['Main Bus B', 11],
            ['Busbar WTG1 LV', 0.69],     
            ['Busbar WTG1 HV', 33], 
            ['Terminal WindPark', 33],
            ['Drilling AC Bus B', 11],
            ['Terminal AC High Drilling B', 11],
            ['Terminal AC Low Drilling B', 0.69],
            ['Terminal DC Drilling B', 1.3],
            ['Drilling DC Bus B', 1.3],
            ['SG Terminal 3', 11],
            ['Terminal ACO_02', 11],
            ['Terminal AC High PEC_VSDc_GEX_02', 11],
            ['Terminal AC Low PEC_VSDc_GEX_02', 3.3],
            ['Terminal DC PEC_VSDc_GEX_02', 6],
            ['Terminal TRA_UTL690_02', 11],
            ['Utility690 Bus B',  0.69],
            ['Terminal TRA_UTL400_02', 0.69],
            ['Utility400 SWBD/BusB', 0.4]

        ],

        'lines': [
            ['name', 'from_bus', 'to_bus', 'length', 'S_n', 'V_n', 'unit', 'R', 'X', 'B'],
            ['Cable WTG1', 'Terminal WindPark', 'Busbar WTG1 HV', 9, 27.95, 33, 'pf', 0.7074, 0.922869, 1.75e-3],
            ['Line GEX_02', 'Main Bus B', 'Terminal AC High PEC_VSDc_GEX_02', 0.25, 9.907, 20, 'pf', 0.01965, 0.02356195, 1.75e-3],
            ['Line ACO_02', 'Main Bus B', 'Terminal ACO_02', 0.25, 6.859, 11, 'pf', 0.038875, 0.02591815, 1.75e-3],
            ['Line Utility B', 'Main Bus B', 'Terminal TRA_UTL690_02', 0.15, 6.859, 11, 'pf', 0.023325, 0.01555089, 1.75e-3],
            ['Line 690 Trafo B', 'Utilty690 Bus B', 'Terminal TRA_UTL400_02', 0.15, 0.6215, 0.69, 'pf', 0.01179, 0.01413717, 1.75e-3],
            ['Line SG3', 'SG Terminal 3', 'Main Bus B', 0.05, 39.63, 11, 'pf', 0.0009825, 0.0039275, 1.75e-3]
            ['Line ACO_02', 'Main Bus B', 'Terminal ACO_02', 0.25, 6.859, 11, 'pf', 0.038875, 0.02591815, 1.75e-3], 
            ['Line Drilling B', 'Main Bus B', 'Drilling AC Bus B', 0.25, 6.859, 11, 'pf', 0.038875, 0.02591815, 1.75e-3],
            ['Line AC Drilling B', 'Drilling AC Bus B', 'Terminal AC High Drilling B', 0.05, 6.859, 11, 'pf', 0.007775, 0.00518363, 1.75e-3],
            ['Line DC Drilling B', 'Terminal DC Drilling B', 'Drilling DC Bus B', 0.01, 6.755, 1.3, 'pf', 0.0001, 0.0002356, 1.75e-3]

        ],

        'transformers': [
            ['name', 'from_bus', 'to_bus', 'S_n', 'V_n_from', 'V_n_to', 'R', 'X'],
            ['Trafo WTG1', 'Busbar WTG1 HV', 'Busbar WTG1 LV', 10, 0.69, 33, 0.5, 0.1],
            ['Trafo WindPark', 'Terminal WindPark', 'Main Bus B', 25, 33, 11, 127, 0.1],
            ['Trafo GEX_02', 'Terminal AC High PEC_VSDc_GEX_02', 'Terminal AC Low PEC_VSDc_GEX_02', 10.5, 11, 3.3, 0, 0],
            ['Terminal AC Low PEC_VSDc_GEX_02', 'Terminal TRA_UTL690_02', 'Utility690 Bus B', 3.3, 11, 0.69, 4, 0.1],
            ['TRA_UTL400_02', 'Utility690 Bus B', 'Utility400 SWBD/BusB', 0.6, 0.69, 0.4, 0.1, 0.1],
            ['TRA_Dril_02', 'Terminal AC High Drilling B', 'Terminal AC Low Drilling B', 3.3, 11, 0.69, 0, 0]

        ],

        'loads': [
            ['name', 'bus', 'P', 'Q', 'model'],
            ['L1', 'B3', 1.05, 0, 'Z'],
            ['L2', 'B3', 1.05, 0, 'Z'],
            ['L3', 'B3', 0.25, 0, 'Z'],
            ['L4', 'B3', 0.25, 0, 'Z']
        ],

        'shunts': [
            ['name', 'bus', 'V_n', 'Q', 'model'],
            ['C1', 'B3', 20, 5, 'Z'],
        ],

        'generators': {
            'GEN': [
                ['name',    'bus',  'S_n',  'V_n',  'P',    'V',    'H',    'D',    'X_d',  'X_q',  'X_d_t',    'X_q_t',    'X_d_st',   'X_q_st',   'T_d0_t',   'T_q0_t',   'T_d0_st',  'T_q0_st'],
                ['G1',      'B1',    28,    11,     15,    1.03,   6.5,    0,      1.8,    1.7,    0.3,        0.55,       0.25,       0.25,       8.0,        0.4,        0.03,       0.05],
                ['G2',      'B1',    28,    11,     15,    1.03,   6.5,    0,      1.8,    1.7,    0.3,        0.55,       0.25,       0.25,       8.0,        0.4,        0.03,       0.05],
                ['G3',      'B1',    28,    11,     15,    1.03,   6.5,    0,      1.8,    1.7,    0.3,        0.55,       0.25,       0.25,       8.0,        0.4,        0.03,       0.05],
                ['G3',      'B1',    28,    11,     15,    1.03,   6.5,    0,      1.8,    1.7,    0.3,        0.55,       0.25,       0.25,       8.0,        0.4,        0.03,       0.05]
            ]
        },

        'vsc': {
            'GridSideConverter_PV': [ 
                ['name',   'bus',    'S_n',      "p_ref_grid",      "v_ref_grid",   'k_p',      'k_v',    'T_p',     'T_v',     'k_pll',   'T_pll',    'T_i',      "i_max"],
                ['WT1',    'B12',      20,         0.0,               1.05,           5,          10,        0.1,        100,        5,        0.1,         0.01,      1.2],
            ]
        },

        'gov': {'TGOV1': [
            ['name',    'gen',  'R',    'D_t',  'V_min',    'V_max',    'T_1',  'T_2',  'T_3'],
            ['GOV1',     'G1',   0.05,   0.02,   0,          1,          0.1,    0.09,   0.2]
            ]
        },

        'avr': {
            'SEXS': [
                ['name',   'gen',      'K',    'T_a',  'T_b',  'T_e',  'E_min',    'E_max'],
                ['AVR1',    'G1',       100,    2.0,    10.0,   0.5,    -3,         3],
            ]
        },

        'pss': {
            'STAB1': [
                ['name',    'gen',  'K',    'T',    'T_1',  'T_2',  'T_3',  'T_4',  'H_lim'],
                ['PSS1',    'G1',   50,     10.0,   0.5,    0.5,    0.05,   0.05,   0.03]
            ]
        }
    }