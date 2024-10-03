def load():
    return {
        'base_mva': 10,
        'f': 50,
        'slack_bus': 'B2',

        'buses': [
            ['name',    'V_n'],
            ['B1',      20],
            ['B2',      20],
        ],

        'lines': [
            ['name',    'from_bus', 'to_bus',   'length',   'S_n',  'V_n',  'unit',     'R',    'X',    'B'],
            ['L1-2',    'B1',       'B2',       25,         20,    20,    'p.u.',     1e-4,   1e-3,   1.75e-3],
        ],

        # 'transformers': [
        #     ['name',    'from_bus', 'to_bus',   'S_n',  'V_n_from', 'V_n_to',   'R',    'X'],
        #     ['T1',      'B1',       'B5',       900,    20,         230,        0,      0.15],
        # ],

        'loads': [
            ['name',    'bus',  'P',    'Q',    'model'],
            ['L1',      'B2',   6,    2,    'Z'],
        ],

        # 'shunts': [
        #     ['name',    'bus',  'V_n',  'Q',    'model'],
        #     ['C1',      'B7',   230,    200,    'Z'],
        #     ['C2',      'B9',   230,    350,    'Z'],
        # ],

        # Changing G1 to parameters found on a PMSG generator
        # https://journals.tubitak.gov.tr/cgi/viewcontent.cgi?article=1759&context=elektrik
        'generators': {
            'GEN': [
                ['name',    'bus',  'S_n',  'V_n',  'P',    'V',    'H',    'D',    'X_d',  'X_q',  'X_d_t',    'X_q_t',    'X_d_st',   'X_q_st',   'T_d0_t',   'T_q0_t',   'T_d0_st',  'T_q0_st'],
                ['G2',      'B2',   15,    20,     7,    1.01,   6.5,    0,      1.8,    1.7,    0.3,        0.55,       0.25,       0.25,       8.0,        0.4,        0.03,       0.05],
            ]
        },

        'gov': {
            'TGOV1': [
                ['name',    'gen',  'R',    'D_t',  'V_min',    'V_max',    'T_1',  'T_2',  'T_3'],
                ['GOV2',     'G2',   0.05,   0.02,   0,          1,          0.1,    0.09,   0.2],
            ]
        },

        'avr': {
            'SEXS': [
                ['name',   'gen',      'K',    'T_a',  'T_b',  'T_e',  'E_min',    'E_max'],
                ['AVR2',    'G2',       100,    2.0,    10.0,   0.5,    -3,         3],
            ]
        },

        'pss': {
            'STAB1': [
                ['name',    'gen',  'K',    'T',    'T_1',  'T_2',  'T_3',  'T_4',  'H_lim'],
                ['PSS2',    'G2',   50,     10.0,   0.5,    0.5,    0.05,   0.05,   0.03],
            ]
        },

        # 'vsc': {
        #     'VSC': [
        # ['name',    'T_pll',    'T_i',  'bus',  'P_K_p',    'P_K_i',    'Q_K_p',    'Q_K_i',    'P_setp',   'Q_setp',   ],
        # ['VSC1',    0.1,        1,      'B8',   0.1,        0.1,        0.1,        0.1,        100,          100],
        #     ]
    
        # }
    }

