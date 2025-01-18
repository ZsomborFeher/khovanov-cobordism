from khovanov import *

def disks_6_1():
    L = Link([(9, 4, 10, 5), (5, 8, 6, 9), (11, 2, 12, 3), (3, 10, 4, 11), (1, 7, 2, 6), (7, 1, 8, 12)])

    S0 = Cobordism(L)
    S0.band_move(-1, (0, 0), (2, 1))
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(-1, (1, 2), (3, 3))
    S1.finish()

    return S0, S1

def disks_8_8():
    L = Link([(12, 16, 13, 15), (16, 12, 17, 11), (14, 6, 15, 5), (4, 14, 5, 13), (6, 4, 7, 3), (2, 8, 3, 7), (1, 10, 2, 11), (9, 0, 10, 1), (17, 8, 0, 9)])

    S0 = Cobordism(L)
    S0.band_move(1, (3, 0), ((2, 1), True), (0, 0))
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(1, (0, 2), ((0, 3), False), (5, 2))
    S1.finish()

    S2 = Cobordism(L)
    S2.band_move(1, (0, 1), ((3, 0), False), (2, 1))
    S2.finish()

    S3 = Cobordism(L)
    S3.band_move(1, (2, 2), ((3, 3), True), (4, 2))
    S3.finish()

    return S0, S1, S2, S3

def disks_8_9():
    L = Link([(16, 8, 17, 7), (13, 6, 14, 7), (5, 12, 6, 13), (10, 2, 11, 1), (0, 10, 1, 9), (8, 0, 9, 17), (11, 4, 12, 5), (15, 2, 16, 3), (3, 14, 4, 15)])

    S0 = Cobordism(L)
    S0.band_move(0, (1, 1), ((8, 2), True), ((7, 1), False), (5, 1))
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(1, (2, 1), ((8, 2), True), (7, 1))
    S1.finish()

    S2 = Cobordism(L)
    S2.band_move(0, (2, 1), ((1, 2), False), ((0, 0), True), (5, 1))
    S2.finish()

    S3 = Cobordism(L)
    S3.band_move(1, (1, 1), ((1, 2), False), (0, 0))
    S3.finish()

    return S0, S1, S2, S3

def disks_9_27():
    L = Link([(12, 15, 13, 16), (9, 19, 10, 18), (17, 9, 18, 8), (16, 5, 17, 6), (2, 11, 3, 12), (6, 1, 7, 2), (19, 11, 0, 10), (0, 7, 1, 8), (14, 4, 15, 3), (4, 14, 5, 13)])

    S0 = Cobordism(L)
    S0.band_move(2, (3, 3), ((3, 0), False), ((3, 1), False), (8, 0))
    S0.finish()
    S0.band_move(0, (9, 0), ((0, 1), False), ((15, 0), False), (4, 0))                  # pickup move
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(2, (5, 3), ((0, 0), True), ((8, 3), True), (8, 0))
    S1.finish()
    S1.band_move(0, (3, 3), ((15, 3), True), ((0, 2), True), (8, 1))                    # pickup move
    S1.finish()

    S2 = Cobordism(L)
    S2.band_move(0, (6, 0), (9, 2))
    S2.finish()
    S2.band_move(1, (6, 1), ((7, 1), True), ((7, 2), True), ((1, 3), True), (1, 0))     # pickup move
    S2.finish()

    S3 = Cobordism(L)
    S3.band_move(0, (1, 0), (4, 2))
    S3.finish()
    S3.band_move(0, (6, 2), ((0, 1), False), ((0, 2), False), (0, 3))                   # pickup move
    S3.finish()

    return S0, S1, S2, S3

def disks_9_41():
    L = Link([(15, 10, 16, 11), (13, 7, 14, 6), (5, 15, 6, 14), (3, 16, 4, 17), (9, 4, 10, 5), (11, 3, 12, 2), (1, 13, 2, 12), (7, 1, 8, 0), (17, 9, 0, 8)])

    S0 = Cobordism(L)
    S0.band_move(1, (1, 3), (5, 0))
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(1, (2, 1), (6, 2))
    S1.finish()

    S2 = Cobordism(L)
    S2.band_move(1, (6, 3), (8, 0))
    S2.finish()

    S3 = Cobordism(L)
    S3.band_move(1, (5, 1), (7, 2))
    S3.finish()

    S4 = Cobordism(L)
    S4.band_move(1, (7, 3), (2, 0))
    S4.finish()

    S5 = Cobordism(L)
    S5.band_move(1, (8, 1), (1, 2))
    S5.finish()

    return S0, S1, S2, S3, S4, S5

def disks_10_3():
    L = Link([(4, 10, 5, 9), (10, 4, 11, 3), (2, 12, 3, 11), (12, 2, 13, 1), (0, 14, 1, 13), (14, 0, 15, 19), (16, 7, 17, 8), (6, 17, 7, 18), (18, 5, 19, 6), (8, 15, 9, 16)])

    S0 = Cobordism(L)
    S0.band_move(1, (0, 3), (4, 0))
    S0.finish()
    
    S1 = Cobordism(L)
    S1.band_move(1, (9, 1), (0, 0))
    S1.finish()

    return S0, S1

def disks_10_123():
    L = Link([(18, 11, 19, 12), (10, 3, 11, 4), (2, 15, 3, 16), (14, 7, 15, 8), (6, 19, 7, 0), (4, 18, 5, 17), (16, 10, 17, 9), (8, 2, 9, 1), (0, 14, 1, 13), (12, 6, 13, 5)])

    S0 = Cobordism(L)
    S0.band_move(1, (7, 2), ((9, 3), False), ((9, 0), False), ((4, 1), True), (2, 2))
    S0.band_move(0, (4, 1), ((2, 2), True), (6, 1))
    S0.finish()

    S1 = Cobordism(L)
    S1.band_move(1, (8, 2), ((5, 3), False), ((5, 0), False), ((0, 1), True), (3, 2))
    S1.band_move(0, (0, 1), ((3, 2), True), (7, 1))
    S1.finish()

    S2 = Cobordism(L)
    S2.band_move(1, (9, 2), ((6, 3), False), ((6, 0), False), ((1, 1), True), (4, 2))
    S2.band_move(0, (1, 1), ((4, 2), True), (8, 1))
    S2.finish()

    S3 = Cobordism(L)
    S3.band_move(1, (5, 2), ((7, 3), False), ((7, 0), False), ((2, 1), True), (0, 2))
    S3.band_move(0, (2, 1), ((0, 2), True), (9, 1))
    S3.finish()

    S4 = Cobordism(L)
    S4.band_move(1, (6, 2), ((8, 3), False), ((8, 0), False), ((3, 1), True), (1, 2))
    S4.band_move(0, (3, 1), ((1, 2), True), (5, 1))
    S4.finish()

    S5 = Cobordism(L)
    S5.band_move(1, (7, 2), ((8, 3), True), ((4, 0), True), ((4, 1), False), (3, 2))
    S5.band_move(0, (7, 0), ((2, 1), False), (0, 2))
    S5.finish()

    S6 = Cobordism(L)
    S6.band_move(1, (8, 2), ((9, 3), True), ((0, 0), True), ((0, 1), False), (4, 2))
    S6.band_move(0, (8, 0), ((3, 1), False), (1, 2))
    S6.finish()

    S7 = Cobordism(L)
    S7.band_move(1, (9, 2), ((5, 3), True), ((1, 0), True), ((1, 1), False), (0, 2))
    S7.band_move(0, (9, 0), ((4, 1), False), (2, 2))
    S7.finish()

    S8 = Cobordism(L)
    S8.band_move(1, (5, 2), ((6, 3), True), ((2, 0), True), ((2, 1), False), (1, 2))
    S8.band_move(0, (5, 0), ((0, 1), False), (3, 2))
    S8.finish()

    S9 = Cobordism(L)
    S9.band_move(1, (6, 2), ((7, 3), True), ((3, 0), True), ((3, 1), False), (2, 2))
    S9.band_move(0, (6, 0), ((1, 1), False), (4, 2))
    S9.finish()

    S10 = Cobordism(L)
    S10.band_move(-1, (8, 2), ((5, 3), True), ((5, 0), False), ((0, 1), False), (4, 2))
    S10.band_move(0, (3, 3), ((8, 2), True), (5, 3))
    S10.finish()

    # S11 - S19 not included as they give the same KJ classes

    return S0, S1, S2, S3, S4, S5, S6, S7, S8, S9, S10


disks = disks_8_9()

printing = False
if printing:
    for S in disks:
        print(S.chi(), S.links[-1].crossings)
        print(S)

printing = True
KJ = [disk.KJ_class(printing) for disk in disks]
print("KJ calculated\n")
M = [disk.mirror().KJ_class(printing) for disk in disks]
print("M calculated\n")

printing = True
new_lines = "\n"
for i, Di in enumerate(disks):
    for j, Dj in enumerate(disks):
        if i < j:
            print(i, j, compare(KJ[i], KJ[j], printing), "KJ", new_lines)
            print(i, j, compare(M[i], M[j], printing), "M", new_lines)
            print("--------------------------", new_lines)
