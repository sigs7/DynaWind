def load():
    return {
        'base_mva': 50,
        'f': 50,
        'slack_bus': 'Grid',

        'buses': [
            ['name',    'V_n'],
            ['Grid', 120],
            ['B120', 120],  # 120 kV bus
            ['B25', 25],  # 25 kV bus
            ['B575T', 25],  # 575 V bus
            ['B575', 0.575],  # 575 V bus
        ],

        'lines': [
            ['name',    'from_bus', 'to_bus',      'length',      'S_n',  'V_n',  'unit',     'R',    'X',    'B'],
            ['Line1-2',  'Grid',       'B120',       10,           100,     120,    'PF',     0.01,   0.4,   0],
            ['Line2',    'B25',       'B575T',       30,           100,     25,     'PF',     0.01,   0.4,   0],
        ],

        'transformers': [
            ['name',    'from_bus', 'to_bus',   'S_n',  'V_n_from',    'V_n_to',      'R',    'X'],
            ['T1',      'B120',       'B25',       50,      120,           25,          0.01,      0.15],
            ['T2',      'B575T',       'B575',      20,     25,         0.575,         0.01,      0.15],
        ],

        'generators': {
            'GEN': [
                ['name',   'bus',  'S_n',  'V_n',    'P',    'V',      'H',    'D',    'X_d',  'X_q',  'X_d_t',    'X_q_t',    'X_d_st',   'X_q_st',   'T_d0_t',   'T_q0_t',   'T_d0_st',  'T_q0_st'],
                ['IB',     'Grid',    1e6,    120,    -5,      1,      1e5,      0,     1.05,   0.66,    0.328,      0.66,       1e-5,      1e-5,         1e5,      10000,          1e5,        1e5],
            ],
        },

        'loads': [
            ['name', 'bus', 'P', 'Q', 'model'],
            ['Load1', 'B120', 15,  0, 'Z'],
        ],

        'vsc': {
            'GridSideConverter': [ 
                ['name',   'bus',    'S_n',      "p_ref_grid",      "q_ref_grid",   'k_p',      'k_q',    'T_p',     'T_q',     'k_pll',   'T_pll',    'T_i',      "i_max"],
                ['WT1',    'B575',      15,         0.6667,               0.0,           5,          1,        0.1,        0.1,        5,        1,         0.01,      1.2],
            ]
        }
    }