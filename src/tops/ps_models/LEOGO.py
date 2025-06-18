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
            ['SG Terminal 3', 11],
            ['Terminal ACO_02', 11],
            ['Terminal AC High PEC_VSDc_GEX_02', 11],
            ['Terminal AC Low PEC_VSDc_GEX_02', 3.3],
            ['Terminal TRA_UTL690_02', 11],
            ['Utility690 Bus B',  0.69],
            ['Terminal TRA_UTL400_02', 0.69],
            ['Utility400 SWBD/BusB', 0.4]

        ],

        'lines': [
            ['name', 'from_bus', 'to_bus', 'length', 'S_n', 'V_n', 'unit', 'R', 'X', 'B'],
            ['Cable WTG1', 'Terminal WindPark', 'Busbar WTG1 HV', 9, 100, 33, 'pf', 0.0786, 0.102541, 0.3],
            ['Line GEX_02', 'Main Bus B', 'Terminal AC High PEC_VSDc_GEX_02', 0.25, 100, 20, 'pf', 0.0786, 0.094247, 0.37],
            ['Line ACO_02', 'Main Bus B', 'Terminal ACO_02', 0.25, 100, 11, 'pf', 0.1555, 0.103672, 0.28],
            ['Line Utility B', 'Main Bus B', 'Terminal TRA_UTL690_02', 0.15, 100, 11, 'pf', 0.1555, 0.103672, 0.28],
            ['Line 690 Trafo B', 'Utilty690 Bus B', 'Terminal TRA_UTL400_02', 100, 0.6215, 0.69, 'pf', 0.0786, 0.094247, 0.37],
            ['Line SG3', 'SG Terminal 3', 'Main Bus B', 0.05, 100, 11, 'pf', 0.0786, 0.094247, 0.37],
            ['Line ACO_02', 'Main Bus B', 'Terminal ACO_02', 0.25, 100, 11, 'pf', 0.1555, 0.103672, 0.28], 
            ['Line Drilling B', 'Main Bus B', 'Drilling AC Bus B', 0.25, 100, 11, 'pf', 0.1555, 0.103672, 0.28],
            ['Line AC Drilling B', 'Drilling AC Bus B', 'Terminal AC High Drilling B', 0.05, 100, 11, 'pf', 0.1555, 0.103672, 0.28]

        ],

        'transformers': [
            ['name', 'from_bus', 'to_bus', 'S_n', 'V_n_from', 'V_n_to', 'R', 'X'],
            ['Trafo WTG1', 'Busbar WTG1 HV', 'Busbar WTG1 LV', 10, 0.69, 33, 0.005, 0.0597913],
            ['Trafo WindPark', 'Terminal WindPark', 'Main Bus B', 25, 33, 11, 0.002, 0.08997778],
            ['Trafo GEX_02', 'Terminal AC High PEC_VSDc_GEX_02', 'Terminal AC Low PEC_VSDc_GEX_02', 10.5, 11, 3.3, 0.0047619, 0.05981074],
            ['TRA_UTL690_02', 'Terminal TRA_UTL690_02', 'Utility690 Bus B', 3.3, 11, 0.69, 0.00606061, 0.1098329],
            ['TRA_UTL400_02', 'Utility690 Bus B', 'Utility400 SWBD/BusB', 0.6, 0.69, 0.4, 0.00833333, 0.05941848],
            ['TRA_Dril_02', 'Terminal AC High Drilling B', 'Terminal AC Low Drilling B', 3.3, 11, 0.69, 0.00606061, 0.1098329]

        ],

        'loads': [
            ['name', 'bus', 'P', 'Q', 'model'],
            ['DC Load - VSD_GEX_02', 'Terminal AC Low PEC_VSDc_GEX_02', 9.7, 0, 'Z'],
            ['LOD2 - Utility Load 690V B', 'Utility690 Bus B', 1.049999, 0, 'Z'],
            ['LOD4 - Utility Load 400V B', 'Utility400 SWBD/BusB', 0.25, 0, 'Z'],
            ['ASM_ACO_02', 'Main Bus B', 0.5, 0.4, 'Z'],
            ['ASM 690V B', 'Utility690 Bus B', 0.2, 0.1, 'Z']
    
        ],

        'shunts': [
            ['name', 'bus', 'V_n', 'Q', 'model'],
            ['LCL-Kondensator Drilling B', 'Terminal AC Low Drilling B', 0.69, 0.139, 'Z'],
            ['LCL-Kondensator GEX_02','Terminal AC Low PEC_VSDc_GEX_02', 3.3, 0.435, 'Z']
        ],

        'generators': {
            'GEN': [
                ['name',    'bus',  'S_n',  'V_n',  'P',    'V',    'H',    'D',    'X_d',  'X_q',  'X_d_t',    'X_q_t',    'X_d_st',   'X_q_st',   'T_d0_t',   'T_q0_t',   'T_d0_st',  'T_q0_st'],
                ['Synchronous Generator 3', 'Main Bus B', 28, 11, 20, 1.02, 6.5, 0, 2.33, 2.1, 0.173, 0, 0.159, 0.159, 0.822, 0, 0.03, 0.013]
            ]
        },

        'gov': {'TGOV1': [ ## GAST -> TGOV1
            ['name',    'gen',  'R',    'D_t',  'V_min',    'V_max',    'T_1',  'T_2',  'T_3'],
            ['GOV1',     'G1',   0.05,   0,   0, 1,    0.4,    3,   0]
            ]
        },

        'avr': {
            'SEXS': [ ## ext_ST
                ['name',   'gen',      'K',    'T_a',  'T_b',  'T_e',  'E_min',    'E_max'],
                ['Static Excitation System', 'Synchronous Generator 3', 50, 0.02, 0, 0.02, -1, 5],
            ]
        }
    }


