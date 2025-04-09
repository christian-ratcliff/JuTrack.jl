# JuTrack lattice file converted from FLAME format
# Original file: BDS_124Xe_3cs_short_v1.lat
# Generated on: 2025-04-09T11:29:39.564
# Mass number: 124.0
# Coordinate system conversion: FLAME [mm, rad, mm, rad, rad, MeV/u] -> JuTrack [m, px, m, py, m, dp]
# Relativistic factors: beta=0.5946004840620179, gamma=1.2437481314969263, beta*gamma=0.7395332410393027
# RF frequency: 8.05e7 Hz, wavelength: 3.724129913043478 m

include("../../JuTrack.jl")
using .JuTrack
using LinearAlgebra
using Random
using Distributions
using CairoMakie
using Statistics

function create_lattice()
    lattice = []

    # Define elements
    LS3_WD06_BPM_D4699 = MARKER(name="LS3_WD06:BPM_D4699")
    push!(lattice, LS3_WD06_BPM_D4699)
    drift_499 = DRIFT(name="drift_499", len=0.999498)
    push!(lattice, drift_499)
    LS3_BTS_DCH_D4709 = ORBTRIM(name="LS3_BTS:DCH_D4709", realpara=true, tm_xkick=-0.0005735004899999999)
    push!(lattice, LS3_BTS_DCH_D4709)
    LS3_BTS_DCV_D4709 = ORBTRIM(name="LS3_BTS:DCV_D4709", realpara=true, tm_ykick=0.000551383048)
    push!(lattice, LS3_BTS_DCV_D4709)
    drift_500 = DRIFT(name="drift_500", len=0.2945)
    push!(lattice, drift_500)
    LS3_BTS_QV_D4713 = KQUAD(name="LS3_BTS:QV_D4713", k1=-7.073181959999999, len=0.261)
    push!(lattice, LS3_BTS_QV_D4713)
    drift_501 = DRIFT(name="drift_501", len=0.239)
    push!(lattice, drift_501)
    LS3_BTS_QH_D4718 = KQUAD(name="LS3_BTS:QH_D4718", k1=7.84223811, len=0.261)
    push!(lattice, LS3_BTS_QH_D4718)
    drift_502 = DRIFT(name="drift_502", len=0.7016340000000001)
    push!(lattice, drift_502)
    LS3_BTS_GV_D4726 = MARKER(name="LS3_BTS:GV_D4726")
    push!(lattice, LS3_BTS_GV_D4726)
    drift_503 = DRIFT(name="drift_503", len=2.421837)
    push!(lattice, drift_503)
    LS3_BTS_DCH_D4750 = ORBTRIM(name="LS3_BTS:DCH_D4750", realpara=true, tm_xkick=-0.00011964992)
    push!(lattice, LS3_BTS_DCH_D4750)
    LS3_BTS_DCV_D4750 = ORBTRIM(name="LS3_BTS:DCV_D4750", realpara=true, tm_ykick=0.00035710039199999996)
    push!(lattice, LS3_BTS_DCV_D4750)
    drift_504 = DRIFT(name="drift_504", len=0.243209)
    push!(lattice, drift_504)
    LS3_BTS_BPM_D4753 = MARKER(name="LS3_BTS:BPM_D4753")
    push!(lattice, LS3_BTS_BPM_D4753)
    drift_505 = DRIFT(name="drift_505", len=0.05232)
    push!(lattice, drift_505)
    LS3_BTS_QV_D4755 = KQUAD(name="LS3_BTS:QV_D4755", k1=-3.9165533400000005, len=0.261)
    push!(lattice, LS3_BTS_QV_D4755)
    drift_506 = DRIFT(name="drift_506", len=0.239)
    push!(lattice, drift_506)
    LS3_BTS_QH_D4760 = KQUAD(name="LS3_BTS:QH_D4760", k1=1.80802314, len=0.261)
    push!(lattice, LS3_BTS_QH_D4760)
    drift_507 = DRIFT(name="drift_507", len=0.826117)
    push!(lattice, drift_507)
    LS3_BTS_BPM_D4769 = MARKER(name="LS3_BTS:BPM_D4769")
    push!(lattice, LS3_BTS_BPM_D4769)
    drift_508 = DRIFT(name="drift_508", len=0.145282)
    push!(lattice, drift_508)
    LS3_BTS_PM_D4771 = MARKER(name="LS3_BTS:PM_D4771")
    push!(lattice, LS3_BTS_PM_D4771)
    drift_509 = DRIFT(name="drift_509", len=0.825101)
    push!(lattice, drift_509)
    LS3_BTS_DCH_D4779 = ORBTRIM(name="LS3_BTS:DCH_D4779", realpara=true, tm_xkick=-0.00010691387)
    push!(lattice, LS3_BTS_DCH_D4779)
    LS3_BTS_DCV_D4779 = ORBTRIM(name="LS3_BTS:DCV_D4779", realpara=true, tm_ykick=-0.0009388861279999999)
    push!(lattice, LS3_BTS_DCV_D4779)
    drift_510 = DRIFT(name="drift_510", len=0.2945)
    push!(lattice, drift_510)
    LS3_BTS_QH_D4783 = KQUAD(name="LS3_BTS:QH_D4783", k1=6.22882116, len=0.261)
    push!(lattice, LS3_BTS_QH_D4783)
    drift_511 = DRIFT(name="drift_511", len=1.199075)
    push!(lattice, drift_511)
    LS3_BTS_PM_D4797 = MARKER(name="LS3_BTS:PM_D4797")
    push!(lattice, LS3_BTS_PM_D4797)
    drift_512 = DRIFT(name="drift_512", len=0.789925)
    push!(lattice, drift_512)
    LS3_BTS_QV_D4806 = KQUAD(name="LS3_BTS:QV_D4806", k1=-6.88877451, len=0.261)
    push!(lattice, LS3_BTS_QV_D4806)
    drift_513 = DRIFT(name="drift_513", len=2.006084)
    push!(lattice, drift_513)
    LS3_BTS_PM_D4827 = MARKER(name="LS3_BTS:PM_D4827")
    push!(lattice, LS3_BTS_PM_D4827)
    drift_514 = DRIFT(name="drift_514", len=1.391397)
    push!(lattice, drift_514)
    LS3_BTS_DCH_D4841 = ORBTRIM(name="LS3_BTS:DCH_D4841", realpara=true, tm_xkick=0.00060231174)
    push!(lattice, LS3_BTS_DCH_D4841)
    LS3_BTS_DCV_D4841 = ORBTRIM(name="LS3_BTS:DCV_D4841", realpara=true, tm_ykick=0.000606018168)
    push!(lattice, LS3_BTS_DCV_D4841)
    drift_515 = DRIFT(name="drift_515", len=0.243199)
    push!(lattice, drift_515)
    LS3_BTS_BPM_D4843 = MARKER(name="LS3_BTS:BPM_D4843")
    push!(lattice, LS3_BTS_BPM_D4843)
    drift_516 = DRIFT(name="drift_516", len=0.05232)
    push!(lattice, drift_516)
    LS3_BTS_QH_D4845 = KQUAD(name="LS3_BTS:QH_D4845", k1=7.06167873, len=0.261)
    push!(lattice, LS3_BTS_QH_D4845)
    drift_517 = DRIFT(name="drift_517", len=1.4945)
    push!(lattice, drift_517)
    LS3_BTS_PM_D4862 = MARKER(name="LS3_BTS:PM_D4862")
    push!(lattice, LS3_BTS_PM_D4862)
    drift_518 = DRIFT(name="drift_518", len=0.4945)
    push!(lattice, drift_518)
    LS3_BTS_QV_D4868 = KQUAD(name="LS3_BTS:QV_D4868", k1=-7.06820118, len=0.261)
    push!(lattice, LS3_BTS_QV_D4868)
    drift_519 = DRIFT(name="drift_519", len=0.158432)
    push!(lattice, drift_519)
    LS3_BTS_GV_D4871 = MARKER(name="LS3_BTS:GV_D4871")
    push!(lattice, LS3_BTS_GV_D4871)
    drift_520 = DRIFT(name="drift_520", len=1.564721)
    push!(lattice, drift_520)
    LS3_BTS_BPM_D4886 = MARKER(name="LS3_BTS:BPM_D4886")
    push!(lattice, LS3_BTS_BPM_D4886)
    drift_521 = DRIFT(name="drift_521", len=1.675347)
    push!(lattice, drift_521)
    LS3_BTS_DCH_D4903 = ORBTRIM(name="LS3_BTS:DCH_D4903", realpara=true, tm_xkick=-0.0006071283899999999)
    push!(lattice, LS3_BTS_DCH_D4903)
    LS3_BTS_DCV_D4903 = ORBTRIM(name="LS3_BTS:DCV_D4903", realpara=true, tm_ykick=-0.0009477944880000001)
    push!(lattice, LS3_BTS_DCV_D4903)
    drift_522 = DRIFT(name="drift_522", len=0.2945)
    push!(lattice, drift_522)
    LS3_BTS_QH_D4907 = KQUAD(name="LS3_BTS:QH_D4907", k1=7.06511784, len=0.261)
    push!(lattice, LS3_BTS_QH_D4907)
    drift_523 = DRIFT(name="drift_523", len=1.9889999999999999)
    push!(lattice, drift_523)
    LS3_BTS_QV_D4930 = KQUAD(name="LS3_BTS:QV_D4930", k1=-7.0700986200000004, len=0.261)
    push!(lattice, LS3_BTS_QV_D4930)
    drift_524 = DRIFT(name="drift_524", len=3.397481)
    push!(lattice, drift_524)
    LS3_BTS_DCH_D4965 = ORBTRIM(name="LS3_BTS:DCH_D4965", realpara=true, tm_xkick=0.0)
    push!(lattice, LS3_BTS_DCH_D4965)
    LS3_BTS_DCV_D4965 = ORBTRIM(name="LS3_BTS:DCV_D4965", realpara=true, tm_ykick=0.0)
    push!(lattice, LS3_BTS_DCV_D4965)
    drift_525 = DRIFT(name="drift_525", len=0.243199)
    push!(lattice, drift_525)
    LS3_BTS_BPM_D4968 = MARKER(name="LS3_BTS:BPM_D4968")
    push!(lattice, LS3_BTS_BPM_D4968)
    drift_526 = DRIFT(name="drift_526", len=0.05232)
    push!(lattice, drift_526)
    LS3_BTS_QH_D4969 = KQUAD(name="LS3_BTS:QH_D4969", k1=7.10555703, len=0.261)
    push!(lattice, LS3_BTS_QH_D4969)
    drift_527 = DRIFT(name="drift_527", len=1.9889999999999999)
    push!(lattice, drift_527)
    LS3_BTS_QV_D4992 = KQUAD(name="LS3_BTS:QV_D4992", k1=-7.11955065, len=0.261)
    push!(lattice, LS3_BTS_QV_D4992)
    drift_528 = DRIFT(name="drift_528", len=1.7231509999999999)
    push!(lattice, drift_528)
    LS3_BTS_BPM_D5010 = MARKER(name="LS3_BTS:BPM_D5010")
    push!(lattice, LS3_BTS_BPM_D5010)
    drift_529 = DRIFT(name="drift_529", len=1.675349)
    push!(lattice, drift_529)
    LS3_BTS_DCH_D5027 = ORBTRIM(name="LS3_BTS:DCH_D5027", realpara=true, tm_xkick=0.00011228606)
    push!(lattice, LS3_BTS_DCH_D5027)
    LS3_BTS_DCV_D5027 = ORBTRIM(name="LS3_BTS:DCV_D5027", realpara=true, tm_ykick=-0.000697141592)
    push!(lattice, LS3_BTS_DCV_D5027)
    drift_530 = DRIFT(name="drift_530", len=0.2945)
    push!(lattice, drift_530)
    LS3_BTS_QH_D5031 = KQUAD(name="LS3_BTS:QH_D5031", k1=7.10140638, len=0.261)
    push!(lattice, LS3_BTS_QH_D5031)
    drift_531 = DRIFT(name="drift_531", len=1.9889999999999999)
    push!(lattice, drift_531)
    LS3_BTS_QV_D5054 = KQUAD(name="LS3_BTS:QV_D5054", k1=-7.1070987, len=0.261)
    push!(lattice, LS3_BTS_QV_D5054)
    drift_532 = DRIFT(name="drift_532", len=3.397481)
    push!(lattice, drift_532)
    LS3_BTS_DCH_D5089 = ORBTRIM(name="LS3_BTS:DCH_D5089", realpara=true, tm_xkick=-0.0002995188)
    push!(lattice, LS3_BTS_DCH_D5089)
    LS3_BTS_DCV_D5089 = ORBTRIM(name="LS3_BTS:DCV_D5089", realpara=true, tm_ykick=0.000304822496)
    push!(lattice, LS3_BTS_DCV_D5089)
    drift_533 = DRIFT(name="drift_533", len=0.243199)
    push!(lattice, drift_533)
    LS3_BTS_BPM_D5092 = MARKER(name="LS3_BTS:BPM_D5092")
    push!(lattice, LS3_BTS_BPM_D5092)
    drift_534 = DRIFT(name="drift_534", len=0.05232)
    push!(lattice, drift_534)
    LS3_BTS_QH_D5093 = KQUAD(name="LS3_BTS:QH_D5093", k1=7.1094705000000005, len=0.261)
    push!(lattice, LS3_BTS_QH_D5093)
    drift_535 = DRIFT(name="drift_535", len=1.9889999999999999)
    push!(lattice, drift_535)
    LS3_BTS_QV_D5116 = KQUAD(name="LS3_BTS:QV_D5116", k1=-7.0999833, len=0.261)
    push!(lattice, LS3_BTS_QV_D5116)
    drift_536 = DRIFT(name="drift_536", len=1.723146854)
    push!(lattice, drift_536)
    LS3_BTS_BPM_D5134 = MARKER(name="LS3_BTS:BPM_D5134")
    push!(lattice, LS3_BTS_BPM_D5134)
    drift_537 = DRIFT(name="drift_537", len=1.675353146)
    push!(lattice, drift_537)
    LS3_BTS_DCH_D5151 = ORBTRIM(name="LS3_BTS:DCH_D5151", realpara=true, tm_xkick=-0.00079380953)
    push!(lattice, LS3_BTS_DCH_D5151)
    LS3_BTS_DCV_D5151 = ORBTRIM(name="LS3_BTS:DCV_D5151", realpara=true, tm_ykick=-0.000319609104)
    push!(lattice, LS3_BTS_DCV_D5151)
    drift_538 = DRIFT(name="drift_538", len=0.2945)
    push!(lattice, drift_538)
    LS3_BTS_QH_D5155 = KQUAD(name="LS3_BTS:QH_D5155", k1=7.10365959, len=0.261)
    push!(lattice, LS3_BTS_QH_D5155)
    drift_539 = DRIFT(name="drift_539", len=1.9889999999999999)
    push!(lattice, drift_539)
    LS3_BTS_QV_D5178 = KQUAD(name="LS3_BTS:QV_D5178", k1=-7.10745447, len=0.261)
    push!(lattice, LS3_BTS_QV_D5178)
    drift_540 = DRIFT(name="drift_540", len=3.3985)
    push!(lattice, drift_540)
    LS3_BTS_DCH_D5213 = ORBTRIM(name="LS3_BTS:DCH_D5213", realpara=true, tm_xkick=0.0)
    push!(lattice, LS3_BTS_DCH_D5213)
    LS3_BTS_DCV_D5213 = ORBTRIM(name="LS3_BTS:DCV_D5213", realpara=true, tm_ykick=0.000645820128)
    push!(lattice, LS3_BTS_DCV_D5213)
    drift_541 = DRIFT(name="drift_541", len=0.242179)
    push!(lattice, drift_541)
    LS3_BTS_BPM_D5216 = MARKER(name="LS3_BTS:BPM_D5216")
    push!(lattice, LS3_BTS_BPM_D5216)
    drift_542 = DRIFT(name="drift_542", len=0.05232099999999999)
    push!(lattice, drift_542)
    LS3_BTS_QH_D5218 = KQUAD(name="LS3_BTS:QH_D5218", k1=7.099627529999999, len=0.261)
    push!(lattice, LS3_BTS_QH_D5218)
    drift_543 = DRIFT(name="drift_543", len=1.9889999999999999)
    push!(lattice, drift_543)
    LS3_BTS_QV_D5240 = KQUAD(name="LS3_BTS:QV_D5240", k1=-7.1041339500000005, len=0.261)
    push!(lattice, LS3_BTS_QV_D5240)
    drift_544 = DRIFT(name="drift_544", len=1.723164854)
    push!(lattice, drift_544)
    LS3_BTS_BPM_D5259 = MARKER(name="LS3_BTS:BPM_D5259")
    push!(lattice, LS3_BTS_BPM_D5259)
    drift_545 = DRIFT(name="drift_545", len=1.675335146)
    push!(lattice, drift_545)
    LS3_BTS_DCH_D5275 = ORBTRIM(name="LS3_BTS:DCH_D5275", realpara=true, tm_xkick=-0.00033805594)
    push!(lattice, LS3_BTS_DCH_D5275)
    LS3_BTS_DCV_D5275 = ORBTRIM(name="LS3_BTS:DCV_D5275", realpara=true, tm_ykick=0.0005861023759999999)
    push!(lattice, LS3_BTS_DCV_D5275)
    drift_546 = DRIFT(name="drift_546", len=0.2945)
    push!(lattice, drift_546)
    LS3_BTS_QH_D5280 = KQUAD(name="LS3_BTS:QH_D5280", k1=7.1107749899999995, len=0.261)
    push!(lattice, LS3_BTS_QH_D5280)
    drift_547 = DRIFT(name="drift_547", len=1.9889999999999999)
    push!(lattice, drift_547)
    LS3_BTS_QV_D5302 = KQUAD(name="LS3_BTS:QV_D5302", k1=-7.09713714, len=0.261)
    push!(lattice, LS3_BTS_QV_D5302)
    drift_548 = DRIFT(name="drift_548", len=3.3985)
    push!(lattice, drift_548)
    LS3_BTS_DCH_D5337 = ORBTRIM(name="LS3_BTS:DCH_D5337", realpara=true, tm_xkick=-0.00063958414)
    push!(lattice, LS3_BTS_DCH_D5337)
    LS3_BTS_DCV_D5337 = ORBTRIM(name="LS3_BTS:DCV_D5337", realpara=true, tm_ykick=-0.0008330268799999999)
    push!(lattice, LS3_BTS_DCV_D5337)
    drift_549 = DRIFT(name="drift_549", len=0.242191)
    push!(lattice, drift_549)
    LS3_BTS_BPM_D5340 = MARKER(name="LS3_BTS:BPM_D5340")
    push!(lattice, LS3_BTS_BPM_D5340)
    drift_550 = DRIFT(name="drift_550", len=0.052308999999999994)
    push!(lattice, drift_550)
    LS3_BTS_QH_D5342 = KQUAD(name="LS3_BTS:QH_D5342", k1=7.10994486, len=0.261)
    push!(lattice, LS3_BTS_QH_D5342)
    drift_551 = DRIFT(name="drift_551", len=1.9889999999999999)
    push!(lattice, drift_551)
    LS3_BTS_QV_D5364 = KQUAD(name="LS3_BTS:QV_D5364", k1=-7.1072172899999995, len=0.261)
    push!(lattice, LS3_BTS_QV_D5364)
    drift_552 = DRIFT(name="drift_552", len=1.583223)
    push!(lattice, drift_552)
    LS3_BTS_BPM_D5381 = MARKER(name="LS3_BTS:BPM_D5381")
    push!(lattice, LS3_BTS_BPM_D5381)
    drift_553 = DRIFT(name="drift_553", len=0.736277)
    push!(lattice, drift_553)
    LS3_BTS_DCH_D5389 = ORBTRIM(name="LS3_BTS:DCH_D5389", realpara=true, tm_xkick=0.00022455636)
    push!(lattice, LS3_BTS_DCH_D5389)
    LS3_BTS_DCV_D5389 = ORBTRIM(name="LS3_BTS:DCV_D5389", realpara=true, tm_ykick=-0.000485012592)
    push!(lattice, LS3_BTS_DCV_D5389)
    drift_554 = DRIFT(name="drift_554", len=0.2945)
    push!(lattice, drift_554)
    LS3_BTS_QH_D5393 = KQUAD(name="LS3_BTS:QH_D5393", k1=7.0869384, len=0.261)
    push!(lattice, LS3_BTS_QH_D5393)
    drift_555 = DRIFT(name="drift_555", len=1.2389999999999999)
    push!(lattice, drift_555)
    LS3_BTS_QV_D5408 = KQUAD(name="LS3_BTS:QV_D5408", k1=-7.09440957, len=0.261)
    push!(lattice, LS3_BTS_QV_D5408)
    drift_556 = DRIFT(name="drift_556", len=1.873481)
    push!(lattice, drift_556)
    LS3_BTS_DCH_D5428 = ORBTRIM(name="LS3_BTS:DCH_D5428", realpara=true, tm_xkick=-0.00118876892)
    push!(lattice, LS3_BTS_DCH_D5428)
    LS3_BTS_DCV_D5428 = ORBTRIM(name="LS3_BTS:DCV_D5428", realpara=true, tm_ykick=0.00026289183999999996)
    push!(lattice, LS3_BTS_DCV_D5428)
    drift_557 = DRIFT(name="drift_557", len=0.243199)
    push!(lattice, drift_557)
    LS3_BTS_BPM_D5430 = MARKER(name="LS3_BTS:BPM_D5430")
    push!(lattice, LS3_BTS_BPM_D5430)
    drift_558 = DRIFT(name="drift_558", len=0.05232)
    push!(lattice, drift_558)
    LS3_BTS_QH_D5432 = KQUAD(name="LS3_BTS:QH_D5432", k1=23.487911465234422, len=0.261)
    push!(lattice, LS3_BTS_QH_D5432)
    drift_559 = DRIFT(name="drift_559", len=0.489)
    push!(lattice, drift_559)
    LS3_BTS_QV_D5440 = KQUAD(name="LS3_BTS:QV_D5440", k1=-17.354942289892875, len=0.261)
    push!(lattice, LS3_BTS_QV_D5440)
    drift_560 = DRIFT(name="drift_560", len=0.371224814)
    push!(lattice, drift_560)
    LS3_BTS_BPM_D5445 = MARKER(name="LS3_BTS:BPM_D5445")
    push!(lattice, LS3_BTS_BPM_D5445)
    drift_561 = DRIFT(name="drift_561", len=2.223275186)
    push!(lattice, drift_561)
    BDS_BTS_DCH_D5467 = ORBTRIM(name="BDS_BTS:DCH_D5467", realpara=true, tm_xkick=0.00045996132)
    push!(lattice, BDS_BTS_DCH_D5467)
    BDS_BTS_DCV_D5467 = ORBTRIM(name="BDS_BTS:DCV_D5467", realpara=true, tm_ykick=0.000516426528)
    push!(lattice, BDS_BTS_DCV_D5467)
    drift_562 = DRIFT(name="drift_562", len=0.2945)
    push!(lattice, drift_562)
    BDS_BTS_QH_D5471 = KQUAD(name="BDS_BTS:QH_D5471", k1=15.703851295404386, len=0.261)
    push!(lattice, BDS_BTS_QH_D5471)
    drift_563 = DRIFT(name="drift_563", len=0.489)
    push!(lattice, drift_563)
    BDS_BTS_QV_D5479 = KQUAD(name="BDS_BTS:QV_D5479", k1=-2.1049294985230427, len=0.261)
    push!(lattice, BDS_BTS_QV_D5479)
    drift_564 = DRIFT(name="drift_564", len=0.32889999999999997)
    push!(lattice, drift_564)
    BDS_BTS_DCH_D5496 = ORBTRIM(name="BDS_BTS:DCH_D5496", realpara=true, tm_xkick=-0.000205109268)
    push!(lattice, BDS_BTS_DCH_D5496)
    BDS_BTS_DCV_D5496 = ORBTRIM(name="BDS_BTS:DCV_D5496", realpara=true, tm_ykick=0.000257928564)
    push!(lattice, BDS_BTS_DCV_D5496)
    drift_565 = DRIFT(name="drift_565", len=0.2921758)
    push!(lattice, drift_565)
    BDS_BTS_BPM_D5499 = MARKER(name="BDS_BTS:BPM_D5499")
    push!(lattice, BDS_BTS_BPM_D5499)
    drift_566 = DRIFT(name="drift_566", len=0.052324199999999994)
    push!(lattice, drift_566)
    BDS_BTS_QV_D5501 = KQUAD(name="BDS_BTS:QV_D5501", k1=-4.804583880544453, len=0.261)
    push!(lattice, BDS_BTS_QV_D5501)
    drift_567 = DRIFT(name="drift_567", len=0.489)
    push!(lattice, drift_567)
    BDS_BTS_QH_D5509 = KQUAD(name="BDS_BTS:QH_D5509", k1=4.793709029952756, len=0.261)
    push!(lattice, BDS_BTS_QH_D5509)
    drift_568 = DRIFT(name="drift_568", len=0.268080814)
    push!(lattice, drift_568)
    BDS_BTS_BPM_D5513 = MARKER(name="BDS_BTS:BPM_D5513")
    push!(lattice, BDS_BTS_BPM_D5513)
    drift_569 = DRIFT(name="drift_569", len=0.145282)
    push!(lattice, drift_569)
    BDS_BTS_PM_D5514 = MARKER(name="BDS_BTS:PM_D5514")
    push!(lattice, BDS_BTS_PM_D5514)
    drift_570 = DRIFT(name="drift_570", len=0.705851186)
    push!(lattice, drift_570)
    BDS_BTS_BCM_D5521 = MARKER(name="BDS_BTS:BCM_D5521")
    push!(lattice, BDS_BTS_BCM_D5521)
    drift_571 = DRIFT(name="drift_571", len=2.571083)
    push!(lattice, drift_571)
    BDS_BTS_GV_D5547 = MARKER(name="BDS_BTS:GV_D5547")
    push!(lattice, BDS_BTS_GV_D5547)
    drift_572 = DRIFT(name="drift_572", len=0.356979)
    push!(lattice, drift_572)
    BDS_BTS_QV_D5552 = KQUAD(name="BDS_BTS:QV_D5552", k1=-15.048493687224616, len=0.261)
    push!(lattice, BDS_BTS_QV_D5552)
    drift_573 = DRIFT(name="drift_573", len=0.489)
    push!(lattice, drift_573)
    BDS_BTS_QH_D5559 = KQUAD(name="BDS_BTS:QH_D5559", k1=19.82095506614822, len=0.261)
    push!(lattice, BDS_BTS_QH_D5559)
    drift_574 = DRIFT(name="drift_574", len=0.2195)
    push!(lattice, drift_574)
    BDS_BTS_DCH_D5563 = ORBTRIM(name="BDS_BTS:DCH_D5563", realpara=true, tm_xkick=-0.0005881375899999999)
    push!(lattice, BDS_BTS_DCH_D5563)
    BDS_BTS_DCV_D5563 = ORBTRIM(name="BDS_BTS:DCV_D5563", realpara=true, tm_ykick=0.0)
    push!(lattice, BDS_BTS_DCV_D5563)
    drift_575 = DRIFT(name="drift_575", len=0.234072)
    push!(lattice, drift_575)
    BDS_BTS_BPM_D5565 = MARKER(name="BDS_BTS:BPM_D5565")
    push!(lattice, BDS_BTS_BPM_D5565)
    drift_576 = DRIFT(name="drift_576", len=0.145282)
    push!(lattice, drift_576)
    BDS_BTS_PM_D5567 = MARKER(name="BDS_BTS:PM_D5567")
    push!(lattice, BDS_BTS_PM_D5567)
    drift_577 = DRIFT(name="drift_577", len=0.295646)
    push!(lattice, drift_577)
    BDS_BBS_DH_D5578_1 = SBEND(name="BDS_BBS:DH_D5578_1", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5578_1)
    BDS_BBS_DH_D5578_2 = SBEND(name="BDS_BBS:DH_D5578_2", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5578_2)
    BDS_BBS_DH_D5578_3 = SBEND(name="BDS_BBS:DH_D5578_3", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5578_3)
    BDS_BBS_DH_D5578_4 = SBEND(name="BDS_BBS:DH_D5578_4", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5578_4)
    drift_578 = DRIFT(name="drift_578", len=0.37951480000000004)
    push!(lattice, drift_578)
    BDS_BBS_S_D5606 = DRIFT(name="BDS_BBS:S_D5606", len=0.2)
    push!(lattice, BDS_BBS_S_D5606)
    drift_579 = DRIFT(name="drift_579", len=0.2)
    push!(lattice, drift_579)
    BDS_BBS_QH_D5611 = KQUAD(name="BDS_BBS:QH_D5611", k1=12.715320875, len=0.4)
    push!(lattice, BDS_BBS_QH_D5611)
    drift_580 = DRIFT(name="drift_580", len=0.3)
    push!(lattice, drift_580)
    BDS_BBS_DCH_D5616 = ORBTRIM(name="BDS_BBS:DCH_D5616", realpara=true, tm_xkick=-3.8544e-6)
    push!(lattice, BDS_BBS_DCH_D5616)
    BDS_BBS_DCV_D5616 = ORBTRIM(name="BDS_BBS:DCV_D5616", realpara=true, tm_ykick=9.636e-7)
    push!(lattice, BDS_BBS_DCV_D5616)
    drift_581 = DRIFT(name="drift_581", len=0.3)
    push!(lattice, drift_581)
    BDS_BBS_QV_D5621 = KQUAD(name="BDS_BBS:QV_D5621", k1=-7.6287745, len=0.4)
    push!(lattice, BDS_BBS_QV_D5621)
    drift_582 = DRIFT(name="drift_582", len=0.230148)
    push!(lattice, drift_582)
    BDS_BBS_BPM_D5625 = MARKER(name="BDS_BBS:BPM_D5625")
    push!(lattice, BDS_BBS_BPM_D5625)
    drift_583 = DRIFT(name="drift_583", len=0.270036)
    push!(lattice, drift_583)
    BDS_BBS_PM_D5628 = MARKER(name="BDS_BBS:PM_D5628")
    push!(lattice, BDS_BBS_PM_D5628)
    drift_584 = DRIFT(name="drift_584", len=0.449816)
    push!(lattice, drift_584)
    BDS_BBS_DH_D5641_0 = SBEND(name="BDS_BBS:DH_D5641_0", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5641_0)
    BDS_BBS_DH_D5641_1 = SBEND(name="BDS_BBS:DH_D5641_1", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5641_1)
    BDS_BBS_DH_D5641_2 = SBEND(name="BDS_BBS:DH_D5641_2", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5641_2)
    BDS_BBS_DH_D5641_3 = SBEND(name="BDS_BBS:DH_D5641_3", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5641_3)
    BDS_BBS_DH_D5641_4 = SBEND(name="BDS_BBS:DH_D5641_4", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5641_4)
    drift_585 = DRIFT(name="drift_585", len=0.370384)
    push!(lattice, drift_585)
    BDS_BBS_BPM_D5653 = MARKER(name="BDS_BBS:BPM_D5653")
    push!(lattice, BDS_BBS_BPM_D5653)
    drift_586 = DRIFT(name="drift_586", len=0.044129)
    push!(lattice, drift_586)
    BDS_BBS_PM_D5653 = MARKER(name="BDS_BBS:PM_D5653")
    push!(lattice, BDS_BBS_PM_D5653)
    drift_587 = DRIFT(name="drift_587", len=0.347358)
    push!(lattice, drift_587)
    BDS_BBS_DCH_D5657 = ORBTRIM(name="BDS_BBS:DCH_D5657", realpara=true, tm_xkick=-6.132e-7)
    push!(lattice, BDS_BBS_DCH_D5657)
    BDS_BBS_DCV_D5657 = ORBTRIM(name="BDS_BBS:DCV_D5657", realpara=true, tm_ykick=1.314e-6)
    push!(lattice, BDS_BBS_DCV_D5657)
    drift_588 = DRIFT(name="drift_588", len=0.338129)
    push!(lattice, drift_588)
    BDS_BBS_DH_D5668_0 = SBEND(name="BDS_BBS:DH_D5668_0", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5668_0)
    BDS_BBS_DH_D5668_1 = SBEND(name="BDS_BBS:DH_D5668_1", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5668_1)
    BDS_BBS_DH_D5668_2 = SBEND(name="BDS_BBS:DH_D5668_2", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5668_2)
    BDS_BBS_DH_D5668_3 = SBEND(name="BDS_BBS:DH_D5668_3", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5668_3)
    BDS_BBS_DH_D5668_4 = SBEND(name="BDS_BBS:DH_D5668_4", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5668_4)
    drift_589 = DRIFT(name="drift_589", len=0.343965)
    push!(lattice, drift_589)
    BDS_BBS_BPM_D5680 = MARKER(name="BDS_BBS:BPM_D5680")
    push!(lattice, BDS_BBS_BPM_D5680)
    drift_590 = DRIFT(name="drift_590", len=0.270037)
    push!(lattice, drift_590)
    BDS_BBS_PM_D5683 = MARKER(name="BDS_BBS:PM_D5683")
    push!(lattice, BDS_BBS_PM_D5683)
    drift_591 = DRIFT(name="drift_591", len=0.335998)
    push!(lattice, drift_591)
    BDS_BBS_QV_D5688 = KQUAD(name="BDS_BBS:QV_D5688", k1=-7.622504125, len=0.4)
    push!(lattice, BDS_BBS_QV_D5688)
    drift_592 = DRIFT(name="drift_592", len=0.3)
    push!(lattice, drift_592)
    BDS_BBS_DCH_D5693 = ORBTRIM(name="BDS_BBS:DCH_D5693", realpara=true, tm_xkick=-9.0228e-6)
    push!(lattice, BDS_BBS_DCH_D5693)
    BDS_BBS_DCV_D5693 = ORBTRIM(name="BDS_BBS:DCV_D5693", realpara=true, tm_ykick=0.0006286176)
    push!(lattice, BDS_BBS_DCV_D5693)
    drift_593 = DRIFT(name="drift_593", len=0.3)
    push!(lattice, drift_593)
    BDS_BBS_QH_D5698 = KQUAD(name="BDS_BBS:QH_D5698", k1=12.706415125, len=0.4)
    push!(lattice, BDS_BBS_QH_D5698)
    drift_594 = DRIFT(name="drift_594", len=0.2)
    push!(lattice, drift_594)
    BDS_BBS_S_D5703 = DRIFT(name="BDS_BBS:S_D5703", len=0.2)
    push!(lattice, BDS_BBS_S_D5703)
    drift_595 = DRIFT(name="drift_595", len=0.37951480000000004)
    push!(lattice, drift_595)
    BDS_BBS_DH_D5731_0 = SBEND(name="BDS_BBS:DH_D5731_0", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5731_0)
    BDS_BBS_DH_D5731_1 = SBEND(name="BDS_BBS:DH_D5731_1", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5731_1)
    BDS_BBS_DH_D5731_2 = SBEND(name="BDS_BBS:DH_D5731_2", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5731_2)
    BDS_BBS_DH_D5731_3 = SBEND(name="BDS_BBS:DH_D5731_3", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5731_3)
    BDS_BBS_DH_D5731_4 = SBEND(name="BDS_BBS:DH_D5731_4", e2=0.0, angle=3.5, e1=0.0, len=0.330050488198)
    push!(lattice, BDS_BBS_DH_D5731_4)
    drift_596 = DRIFT(name="drift_596", len=0.255853)
    push!(lattice, drift_596)
    BDS_FFS_BPM_D5742 = MARKER(name="BDS_FFS:BPM_D5742")
    push!(lattice, BDS_FFS_BPM_D5742)
    drift_597 = DRIFT(name="drift_597", len=0.145282)
    push!(lattice, drift_597)
    BDS_FFS_PM_D5743 = MARKER(name="BDS_FFS:PM_D5743")
    push!(lattice, BDS_FFS_PM_D5743)
    drift_598 = DRIFT(name="drift_598", len=0.299265)
    push!(lattice, drift_598)
    BDS_FFS_DCH_D5746 = ORBTRIM(name="BDS_FFS:DCH_D5746", realpara=true, tm_xkick=-0.0025638781699999997)
    push!(lattice, BDS_FFS_DCH_D5746)
    BDS_FFS_DCV_D5746 = ORBTRIM(name="BDS_FFS:DCV_D5746", realpara=true, tm_ykick=0.0)
    push!(lattice, BDS_FFS_DCV_D5746)
    drift_599 = DRIFT(name="drift_599", len=0.1941)
    push!(lattice, drift_599)
    BDS_FFS_QH_D5750 = KQUAD(name="BDS_FFS:QH_D5750", k1=11.51497041, len=0.261)
    push!(lattice, BDS_FFS_QH_D5750)
    drift_600 = DRIFT(name="drift_600", len=0.489)
    push!(lattice, drift_600)
    BDS_FFS_QV_D5757 = KQUAD(name="BDS_FFS:QV_D5757", k1=-8.6416533, len=0.261)
    push!(lattice, BDS_FFS_QV_D5757)
    drift_601 = DRIFT(name="drift_601", len=0.27309419999999995)
    push!(lattice, drift_601)
    BDS_FFS_BPM_D5772 = MARKER(name="BDS_FFS:BPM_D5772")
    push!(lattice, BDS_FFS_BPM_D5772)
    drift_602 = DRIFT(name="drift_602", len=0.145282)
    push!(lattice, drift_602)
    BDS_FFS_PM_D5774 = MARKER(name="BDS_FFS:PM_D5774")
    push!(lattice, BDS_FFS_PM_D5774)
    drift_603 = DRIFT(name="drift_603", len=0.206184)
    push!(lattice, drift_603)
    BDS_FFS_QH_D5777 = KQUAD(name="BDS_FFS:QH_D5777", k1=7.48468926, len=0.261)
    push!(lattice, BDS_FFS_QH_D5777)
    drift_604 = DRIFT(name="drift_604", len=0.2695)
    push!(lattice, drift_604)
    BDS_FFS_DCH_D5781 = ORBTRIM(name="BDS_FFS:DCH_D5781", realpara=true, tm_xkick=0.0)
    push!(lattice, BDS_FFS_DCH_D5781)
    BDS_FFS_DCV_D5781 = ORBTRIM(name="BDS_FFS:DCV_D5781", realpara=true, tm_ykick=0.0)
    push!(lattice, BDS_FFS_DCV_D5781)
    drift_605 = DRIFT(name="drift_605", len=0.2195)
    push!(lattice, drift_605)
    BDS_FFS_QV_D5784 = KQUAD(name="BDS_FFS:QV_D5784", k1=-8.503614540000001, len=0.261)
    push!(lattice, BDS_FFS_QV_D5784)
    drift_606 = DRIFT(name="drift_606", len=0.180737)
    push!(lattice, drift_606)
    BDS_FFS_GV_D5788 = MARKER(name="BDS_FFS:GV_D5788")
    push!(lattice, BDS_FFS_GV_D5788)
    drift_607 = DRIFT(name="drift_607", len=0.123632)
    push!(lattice, drift_607)
    BDS_FFS_BCM_D5789 = MARKER(name="BDS_FFS:BCM_D5789")
    push!(lattice, BDS_FFS_BCM_D5789)
    drift_608 = DRIFT(name="drift_608", len=0.123499)
    push!(lattice, drift_608)
    BDS_FFS_BPM_D5790 = MARKER(name="BDS_FFS:BPM_D5790")
    push!(lattice, BDS_FFS_BPM_D5790)
    drift_609 = DRIFT(name="drift_609", len=0.145282)
    push!(lattice, drift_609)
    BDS_FFS_PM_D5792 = MARKER(name="BDS_FFS:PM_D5792")
    push!(lattice, BDS_FFS_PM_D5792)
    drift_610 = DRIFT(name="drift_610", len=0.23193000000000003)
    push!(lattice, drift_610)
    BDS_FFS_BPM_D5803 = MARKER(name="BDS_FFS:BPM_D5803")
    push!(lattice, BDS_FFS_BPM_D5803)
    drift_611 = DRIFT(name="drift_611", len=0.2367)
    push!(lattice, drift_611)
    BDS_FFS_DCH_D5805 = ORBTRIM(name="BDS_FFS:DCH_D5805", realpara=true, tm_xkick=0.0017569932)
    push!(lattice, BDS_FFS_DCH_D5805)
    BDS_FFS_DCV_D5805 = ORBTRIM(name="BDS_FFS:DCV_D5805", realpara=true, tm_ykick=-0.0017488464)
    push!(lattice, BDS_FFS_DCV_D5805)
    drift_612 = DRIFT(name="drift_612", len=0.275063)
    push!(lattice, drift_612)
    BDS_FFS_QH_D5810 = KQUAD(name="BDS_FFS:QH_D5810", k1=20.273760000000003, len=0.26)
    push!(lattice, BDS_FFS_QH_D5810)
    drift_613 = DRIFT(name="drift_613", len=0.24)
    push!(lattice, drift_613)
    BDS_FFS_QV_D5815 = KQUAD(name="BDS_FFS:QV_D5815", k1=-21.3058104, len=0.4)
    push!(lattice, BDS_FFS_QV_D5815)
    drift_614 = DRIFT(name="drift_614", len=0.12288)
    push!(lattice, drift_614)
    BDS_FFS_BPM_D5818 = MARKER(name="BDS_FFS:BPM_D5818")
    push!(lattice, BDS_FFS_BPM_D5818)
    drift_615 = DRIFT(name="drift_615", len=0.11712)
    push!(lattice, drift_615)
    BDS_FFS_QH_D5821 = KQUAD(name="BDS_FFS:QH_D5821", k1=17.316936000000002, len=0.26)
    push!(lattice, BDS_FFS_QH_D5821)
    drift_616 = DRIFT(name="drift_616", len=0.116366)
    push!(lattice, drift_616)
    BDS_FFS_CLLM_D5824 = DRIFT(name="BDS_FFS:CLLM_D5824", len=0.0738636)
    push!(lattice, BDS_FFS_CLLM_D5824)
    drift_617 = DRIFT(name="drift_617", len=0.1103716)
    push!(lattice, drift_617)
    mirror_shield_2 = DRIFT(name="mirror_shield_2", len=0.0127)
    push!(lattice, mirror_shield_2)
    drift_618 = DRIFT(name="drift_618", len=0.045896)
    push!(lattice, drift_618)
    mirror = DRIFT(name="mirror", len=0.011)
    push!(lattice, mirror)
    drift_619 = DRIFT(name="drift_619", len=0.386459)
    push!(lattice, drift_619)
    Target_Center = MARKER(name="Target_Center")
    push!(lattice, Target_Center)
    USE = MARKER(name="USE")
    push!(lattice, USE)

    # Define FLAME raw matrices
    S0_raw = [
        1.9401120068506905 6.916101751370262e-6 -0.6452826622660077 1.937737875971762e-5 2.5016619575190483e-7 -7.22289629767321e-7 ;
        6.916101751370465e-6 1.227632215382666e-8 7.966019841914919e-5 1.2214490136665444e-8 1.937446677836609e-11 1.4392049726901372e-10 ;
        -0.6452826622660028 7.966019841914866e-5 1.6772745702339231 3.786049860001454e-5 -9.291570665087502e-7 -2.371472911673177e-6 ;
        1.9377378759717306e-5 1.221449013666544e-8 3.786049860001564e-5 1.3727527372393609e-8 8.669554527079634e-11 -7.809804442718887e-10 ;
        2.501661957483989e-7 1.9374466777649283e-11 -9.291570665010154e-7 8.669554527099378e-11 5.397896688007718e-5 4.2791987513052936e-5 ;
        -7.222896297913816e-7 1.4392049726288899e-10 -2.371472911608706e-6 -7.809804442646337e-10 4.279198751305285e-5 0.0003339580248268087 ;
    ]

    S1_raw = [
        1.4591505220806416 4.80759780722691e-6 -0.25740080488810885 -1.5753048728498944e-5 2.5016619575190346e-7 -7.222896297673236e-7 ;
        4.807597807226825e-6 1.113553929141832e-8 6.517280550397521e-5 -9.059142689309891e-9 1.9374466778366148e-11 1.4392049726901385e-10 ;
        -0.2574008048881045 6.51728055039751e-5 2.3112947543955262 4.479549584864729e-5 -9.291570665087705e-7 -2.3714729116731894e-6 ;
        -1.5753048728498907e-5 -9.059142689309881e-9 4.479549584864725e-5 1.2334170148319723e-8 8.669554527079657e-11 -7.809804442718882e-10 ;
        2.5016619574839727e-7 1.937446677764928e-11 -9.291570665010305e-7 8.669554527099342e-11 5.397896688007715e-5 4.279198751305294e-5 ;
        -7.222896297913825e-7 1.43920497262889e-10 -2.3714729116086933e-6 -7.80980444264636e-10 4.279198751305286e-5 0.0003339580248268087 ;
    ]

    S2_raw = [
        8.610934701851642 -0.00029259969095107527 -1.0426799474860453 0.000396584305304612 2.501661957519043e-7 -7.222896297673205e-7 ;
        -0.00029259969095107716 2.0119976600609268e-8 9.453461098576867e-5 -3.284566499545562e-9 1.9374466778366112e-11 1.439204972690139e-10 ;
        -1.0426799474860375 9.453461098576793e-5 1.2714676786183123 1.2290925637784062e-5 -9.291570665087563e-7 -2.3714729116731957e-6 ;
        0.0003965843053046115 -3.2845664995454693e-9 1.2290925637784374e-5 2.8471963564199213e-8 8.66955452707956e-11 -7.809804442718881e-10 ;
        2.501661957483968e-7 1.9374466777649218e-11 -9.291570665010196e-7 8.66955452709927e-11 5.3978966880077166e-5 4.279198751305295e-5 ;
        -7.222896297913841e-7 1.439204972628889e-10 -2.371472911608734e-6 -7.809804442646341e-10 4.279198751305287e-5 0.0003339580248268087 ;
    ]

    # Define RF and beam parameters
    rf_freq = 8.05e7  # RF frequency in Hz
    rf_k = 2 * pi * rf_freq / 299792458.0  # [1/m]
    IonEs = 9.3149432e8  # Nucleon mass [eV/u]
    IonEk = 2.2705e8  # Kinetic energy [eV/u]
    gamma = 1.0 + IonEk / IonEs  # Lorentz factor
    beta = sqrt(1.0 - 1.0 / (gamma * gamma))  # Velocity factor
    
    # Define scaling factors for coordinate transformation
    scaling = [1e-3, 1.0, 1e-3, 1.0, beta/rf_k, 1.0e6/beta/beta/(IonEs+IonEk)]
    
    # Transform FLAME matrices to JuTrack matrices
    # Transform S0
    S0_matrix = diagm(scaling) * ((S0_raw + S0_raw') ./ 2) * diagm(scaling)'  # Apply scaling and symmetrize

    # Transform S1
    S1_matrix = diagm(scaling) * ((S1_raw + S1_raw') ./ 2) * diagm(scaling)'  # Apply scaling and symmetrize

    # Transform S2
    S2_matrix = diagm(scaling) * ((S2_raw + S2_raw') ./ 2) * diagm(scaling)'  # Apply scaling and symmetrize

    # Create beam objects for each charge state
    # Generate beam1 using eigendecomposition
    Random.seed!(42)  # Ensure reproducibility
    
    # Get eigendecomposition of the transformed matrix
    lam1, u1 = eigen(S0_matrix)
    
    if any(lam1 .<= 0)
        println("Warning: Detected negative eigenvalues in the matrix")
        # Set threshold relative to the largest eigenvalue
        min_eigenvalue = max(1e-20, maximum(lam1) * 1e-6)  
        # Clip negative eigenvalues
        lam1_clipped = max.(lam1, min_eigenvalue)  
        # Reconstruct matrix with clipped eigenvalues but same eigenvectors
        S_reg = u1 * Diagonal(lam1_clipped) * u1'
        lam1, u1 = eigen(S_reg)
    end

    # Generate random particles 
    nparticles = 10000000
    dis1 = Matrix{Float64}(undef, nparticles, 6)
    for d in 1:6
        dis1[:, d] .= randn(nparticles)
    end
    
    # Calculate the initial 2nd moment matrix 
    moment2nd1 = zeros(6, 6)
    for d1 in 1:6
        for d2 in 1:6
            moment2nd1[d1, d2] = mean(dis1[:, d1] .* dis1[:, d2])
        end
    end
    
    # Create transformation matrix from eigendecomposition
    transformation1 = u1 * Diagonal(sqrt.(lam1)) * u1'
    
    # Apply transformation to get desired covariance
    dis1 = dis1 * transformation1'
    
    # Calculate the final 2nd moment matrix to verify transformation
    for d1 in 1:6
        for d2 in 1:6
            moment2nd1[d1, d2] = mean(dis1[:, d1] .* dis1[:, d2])
        end
    end
    
    # Calculate and display the ratio to verify accuracy
    ratio1 = moment2nd1 ./ S0_matrix
    println("Beam 1 second moment ratio:")
    display(ratio1)
    
    # Add centroid
    dis1 .+= reshape([0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1, 6)
    
    # Create beam with the generated particles - using a subset for efficiency
    beam1 = Beam(r=dis1[1:10000,:], np=10000, energy=2.81542e10, charge=50.0, mass=1.1550529568e11)
    get_centroid!(beam1)
    get_emittance!(beam1)

    # Generate beam2 using eigendecomposition
    Random.seed!(43)  # Ensure reproducibility
    
    # Get eigendecomposition of the transformed matrix
    lam2, u2 = eigen(S1_matrix)
    
    if any(lam2 .<= 0)
        println("Warning: Detected negative eigenvalues in the matrix")
        # Set threshold relative to the largest eigenvalue
        min_eigenvalue = max(1e-20, maximum(lam2) * 1e-6)  
        # Clip negative eigenvalues
        lam1_clipped = max.(lam2, min_eigenvalue)  
        # Reconstruct matrix with clipped eigenvalues but same eigenvectors
        S_reg = u2 * Diagonal(lam1_clipped) * u2'
        lam2, u2 = eigen(S_reg)
    end

    # Generate random particles 
    nparticles = 10000000
    dis2 = Matrix{Float64}(undef, nparticles, 6)
    for d in 1:6
        dis2[:, d] .= randn(nparticles)
    end
    
    # Calculate the initial 2nd moment matrix 
    moment2nd2 = zeros(6, 6)
    for d1 in 1:6
        for d2 in 1:6
            moment2nd2[d1, d2] = mean(dis2[:, d1] .* dis2[:, d2])
        end
    end
    
    # Create transformation matrix from eigendecomposition
    transformation2 = u2 * Diagonal(sqrt.(lam2)) * u2'
    
    # Apply transformation to get desired covariance
    dis2 = dis2 * transformation2'
    
    # Calculate the final 2nd moment matrix to verify transformation
    for d1 in 1:6
        for d2 in 1:6
            moment2nd2[d1, d2] = mean(dis2[:, d1] .* dis2[:, d2])
        end
    end
    
    # Calculate and display the ratio to verify accuracy
    ratio2 = moment2nd2 ./ S1_matrix
    println("Beam 2 second moment ratio:")
    display(ratio2)
    
    # Add centroid
    dis2 .+= reshape([0.0, 0.0, 0.0, 0.0, 0.0, 0.0007812442053822318], 1, 6)
    
    # Create beam with the generated particles - using a subset for efficiency
    beam2 = Beam(r=dis2[1:10000,:], np=10000, energy=2.81542e10, charge=49.0, mass=1.1550529568e11)
    get_centroid!(beam2)
    get_emittance!(beam2)

    # Generate beam3 using eigendecomposition
    Random.seed!(44)  # Ensure reproducibility
    
    # Get eigendecomposition of the transformed matrix
    lam3, u3 = eigen(S2_matrix)
    
    if any(lam3 .<= 0)
        println("Warning: Detected negative eigenvalues in the matrix")
        # Set threshold relative to the largest eigenvalue
        min_eigenvalue = max(1e-20, maximum(lam3) * 1e-6)  
        # Clip negative eigenvalues
        lam1_clipped = max.(lam3, min_eigenvalue)  
        # Reconstruct matrix with clipped eigenvalues but same eigenvectors
        S_reg = u3 * Diagonal(lam1_clipped) * u3'
        lam3, u3 = eigen(S_reg)
    end

    # Generate random particles 
    nparticles = 10000000
    dis3 = Matrix{Float64}(undef, nparticles, 6)
    for d in 1:6
        dis3[:, d] .= randn(nparticles)
    end
    
    # Calculate the initial 2nd moment matrix 
    moment2nd3 = zeros(6, 6)
    for d1 in 1:6
        for d2 in 1:6
            moment2nd3[d1, d2] = mean(dis3[:, d1] .* dis3[:, d2])
        end
    end
    
    # Create transformation matrix from eigendecomposition
    transformation3 = u3 * Diagonal(sqrt.(lam3)) * u3'
    
    # Apply transformation to get desired covariance
    dis3 = dis3 * transformation3'
    
    # Calculate the final 2nd moment matrix to verify transformation
    for d1 in 1:6
        for d2 in 1:6
            moment2nd3[d1, d2] = mean(dis3[:, d1] .* dis3[:, d2])
        end
    end
    
    # Calculate and display the ratio to verify accuracy
    ratio3 = moment2nd3 ./ S2_matrix
    println("Beam 3 second moment ratio:")
    display(ratio3)
    
    # Add centroid
    dis3 .+= reshape([0.0, 0.0, 0.0, 0.0, 0.0, -0.0012939357151643517], 1, 6)
    
    # Create beam with the generated particles - using a subset for efficiency
    beam3 = Beam(r=dis3[1:10000,:], np=10000, energy=2.81542e10, charge=51.0, mass=1.1550529568e11)
    get_centroid!(beam3)
    get_emittance!(beam3)

    return lattice, beam1, beam2, beam3
end

# Run the simulation

function propagate_beam(lattice, beam, np)
    # Create expanded lattice
    expanded_lattice = []
    expanded_indices = []  # Track which original element each expanded element belongs to
    
    for (i, element) in enumerate(lattice)
        if element.len > 0
            # Calculate segments needed
            segments_by_count = 10
            segments_by_length = ceil(Int, element.len / 0.1)
            num_segments = max(segments_by_count, segments_by_length)
            segment_length = element.len / num_segments
            
            # Create smaller elements
            for j in 1:num_segments
                small_element = deepcopy(element)
                small_element.len = segment_length
                push!(expanded_lattice, small_element)
                push!(expanded_indices, i)  # Record original index
            end
        else
            # Zero-length elements
            push!(expanded_lattice, element)
            push!(expanded_indices, i)
        end
    end
    
    # Propagate through expanded lattice
    n_expanded = length(expanded_lattice)
    expanded_beam_rms = zeros(n_expanded, 6)
    expanded_floor_distance = zeros(n_expanded)
    expanded_twi = zeros(n_expanded, 9)
    flat_particles = zeros(6 * np)
    
    # Propagate through expanded lattice
    for i in eachindex(expanded_lattice)
        expanded_twi[i, :] .= twiss_beam(beam)
        flat_particles .= collect(Iterators.flatten(eachrow(beam.r)))
        
        # Pass particles through element
        pass!(expanded_lattice[i], flat_particles, np, beam)
        beam.r = reshape(flat_particles, 6, np)'
        
        # Store RMS values
        for dim in 1:6
            expanded_beam_rms[i, dim] = sqrt(mean(beam.r[:,dim].^2))
        end
        
        # Calculate floor distance
        if i == 1
            expanded_floor_distance[i] = expanded_lattice[i].len
        else
            expanded_floor_distance[i] = expanded_floor_distance[i-1] + expanded_lattice[i].len
        end
    end
    
    # Return everything needed
    return expanded_beam_rms, expanded_floor_distance, beam, expanded_twi, 
        lattice, expanded_lattice, expanded_indices
end

function run_simulation()
    lattice, beam1, beam2, beam3 = create_lattice()
    
    beam1_rms, floor_distance, beam1, twi1, original_lattice, _, expanded_indices1 = 
        propagate_beam(lattice, beam1, beam1.np)
    
    beam2_rms, _, beam2, twi2, _, _, expanded_indices2 = 
        propagate_beam(lattice, beam2, beam2.np)
    
    beam3_rms, _, beam3, twi3, _, _, expanded_indices3 = 
        propagate_beam(lattice, beam3, beam3.np)

    return original_lattice, beam1, beam2, beam3, 
        beam1_rms, beam2_rms, beam3_rms, 
        twi1, twi2, twi3, floor_distance, expanded_indices1
end

# plotting wrapper function
function plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, 
                                    twi1, twi2, twi3, floor_distance,
                                    original_lattice, expanded_indices, 
                                    combined=false)
    # Filter expanded data up to the specified original element
    mask = expanded_indices .<= end_ele
    
    # Use the filtered data for plots but original lattice for floor layout
    return plot_multibeam_data(
        floor_distance[mask], 
        beam1_rms[mask,:], 
        beam2_rms[mask,:], 
        beam3_rms[mask,:], 
        twi1[mask,:], 
        twi2[mask,:], 
        twi3[mask,:], 
        original_lattice[1:end_ele],  # Use original lattice for floor plot
        combined=combined
    )
end

function visualize_beam_properties(beam::Beam)
    # Create a multi-panel figure
    fig = Figure(size=(900, 600))
    
    # Phase space plots (x-px, y-py, z-dp)
    ax1 = Axis(fig[1, 1], title="x-px phase space", xlabel="x [mm]", ylabel="px")
    # Convert m to mm for x-axis
    scatter!(ax1, beam.r[:,1] .* 1000, beam.r[:,2], markersize=1, strokewidth=0, alpha=0.5)
    
    ax2 = Axis(fig[1, 2], title="y-py phase space", xlabel="y [mm]", ylabel="py")
    # Convert m to mm for y-axis
    scatter!(ax2, beam.r[:,3] .* 1000, beam.r[:,4], markersize=1, strokewidth=0, alpha=0.5)
    
    ax3 = Axis(fig[1, 3], title="z-dp phase space", xlabel="z [mm]", ylabel="dp/p")
    # Convert m to mm for z-axis
    scatter!(ax3, beam.r[:,5] .* 1000, beam.r[:,6], markersize=1, strokewidth=0, alpha=0.5)
    
    # Projections (histograms)
    ax4 = Axis(fig[2, 1], title="x distribution", xlabel="x [mm]")
    # Convert m to mm for x-axis
    hist!(ax4, beam.r[:,1] .* 1000, bins=50, color=(:blue, 0.7))
    
    ax5 = Axis(fig[2, 2], title="y distribution", xlabel="y [mm]")
    # Convert m to mm for y-axis
    hist!(ax5, beam.r[:,3] .* 1000, bins=50, color=(:blue, 0.7))
    
    ax6 = Axis(fig[2, 3], title="z distribution", xlabel="z [mm]")
    # Convert m to mm for z-axis
    hist!(ax6, beam.r[:,5] .* 1000, bins=50, color=(:blue, 0.7))
    
    # Return the figure
    return fig
end

function plot_multibeam_data(floor_length, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, lattice; combined = true)
    if combined != true    
        # Define beam colors and labels (shared across all plots)
        beam_colors = [:blue, :orange, :green]
        beam_labels = [L"^{124}Xe^{49+}", L"^{124}Xe^{50+}", L"^{124}Xe^{51+}"]
        
        # Define element colors and heights (shared across all plots)
        element_colors = Dict(
            "SBEND" => :purple,
            "KQUAD" => :red,
            "ORBTRIM" => :green,
            "MARKER" => :gray,
            "DRIFT" => :white,
            "default" => :gray
        )
        
        element_heights = Dict(
            "SBEND" => 0.8,
            "KQUAD" => 0.7,
            "ORBTRIM" => 1.0,
            "MARKER" => 0.4,
            "DRIFT" => 0.3,
            "default" => 0.5
        )
        
        # Helper function to create a standard plot structure
        function create_plot(top_data1, top_data2, top_data3, bottom_data1, bottom_data2, bottom_data3, 
                            top_label, bottom_label, title)
            
            f = Figure(size = (1200, 400))
            gl = f[1, 1] = GridLayout()
            legend_layout = f[1, 2] = GridLayout()
            
            # Top plot
            ax_top = Axis(gl[1, 1], xlabel = "", ylabel = top_label)
            lines!(ax_top, floor_length, top_data1, color = beam_colors[1], linewidth = 2)
            lines!(ax_top, floor_length, top_data2, color = beam_colors[2], linewidth = 2)
            lines!(ax_top, floor_length, top_data3, color = beam_colors[3], linewidth = 2)
            hidexdecorations!(ax_top)
            
            # Create floorline plot in the middle
            floor_ax = Axis(gl[2, 1], 
                        xlabel = "",
                        ylabel = "",
                        yticklabelsvisible = false,
                        yticksvisible = false)
            
            # Track element types for legend
            element_types = Dict{String, Any}()
            
            # Draw the floorline
            curr_pos = 0.0
            for ele in lattice
                # Get element type
                ele_type = string(typeof(ele).name.name)
                
                # Determine color and height
                color = get(element_colors, ele_type, element_colors["default"])
                height = get(element_heights, ele_type, element_heights["default"])
                
                # Make sure elements don't overlap by using a minimum width
                ele_width = ele.len #max(ele.len, 0.05)
                
                # Draw rectangle for element based on element type and properties
                if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
                    # Use the original quad height
                    quad_height = element_heights["KQUAD"]
                    drift_height = element_heights["DRIFT"]
                    
                    # Calculate offset to align with other elements
                    # This will make it overlap the x-axis slightly to maintain alignment
                    if ele.k1 > 0
                        # For focusing quads: align bottom with the bottom of drift elements
                        y_pos = -drift_height/2
                    else
                        # For defocusing quads: align top with the top of drift elements
                        y_pos = drift_height/2 - quad_height
                    end
                    
                    rect = Rect(curr_pos, y_pos, ele_width, quad_height)
                elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
                    # Use the original sextupole height
                    sext_height = get(element_heights, "KSEXT", element_heights["default"])
                    drift_height = element_heights["DRIFT"]
                    
                    # Calculate offset to align with other elements
                    if ele.k2 > 0
                        # For focusing sextupoles: align bottom with the bottom of drift elements
                        y_pos = -drift_height/2
                    else
                        # For defocusing sextupoles: align top with the top of drift elements
                        y_pos = drift_height/2 - sext_height
                    end
                    
                    rect = Rect(curr_pos, y_pos, ele_width, sext_height)
                else
                    # All other elements centered on x-axis as before
                    rect = Rect(curr_pos, -height/2, ele_width, height)
                end
                
                element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
                
                # Track for legend (only if not already tracked)
                if !haskey(element_types, ele_type)
                    element_types[ele_type] = element
                end
                
                # Update position - ensure we advance by at least the element width
                # to avoid overlapping due to position calculations
                curr_pos += ele_width
            end
            
            hidexdecorations!(floor_ax)
            
            # Bottom plot
            ax_bottom = Axis(gl[3, 1], xlabel = "z [m]", ylabel = bottom_label, yreversed = true)
            lines!(ax_bottom, floor_length, bottom_data1, color = beam_colors[1], linewidth = 2)
            lines!(ax_bottom, floor_length, bottom_data2, color = beam_colors[2], linewidth = 2)
            lines!(ax_bottom, floor_length, bottom_data3, color = beam_colors[3], linewidth = 2)
            
            # Link x axes
            linkxaxes!(ax_top, floor_ax, ax_bottom)
            
            # Add beam legend
            legend_entries = [
                LineElement(color = c, linewidth = 2) for c in beam_colors
            ]
            
            leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
            
            # Add element legend (filtering out MARKER and DRIFT)
            element_entries = []
            element_names = []
            sorted_names = sort(collect(keys(element_types)))
            for name in sorted_names
                if !(name in ["MARKER", "DRIFT"])
                    element = element_types[name]
                    push!(element_entries, element)
                    push!(element_names, name)
                end
            end
            
            if !isempty(element_entries)
                element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
            end
            
            # Adjust spacing between legends
            rowgap!(legend_layout, 5)
            
            # Set row sizes
            rowsize!(gl, 1, 100)  # Top plot
            rowsize!(gl, 2, 30)   # Elements plot (smaller)
            rowsize!(gl, 3, 100)  # Bottom plot
            
            # Set title
            # f.title = title
            
            return f
        end
        
        # Create the four plots
        rms_plot = create_plot(
            beam1_rms[:,1] .* 1e3,
            beam2_rms[:,1] .* 1e3,
            beam3_rms[:,1] .* 1e3,
            beam1_rms[:,3] .* 1e3,
            beam2_rms[:,3] .* 1e3,
            beam3_rms[:,3] .* 1e3,
            "RMS X [mm]",
            "RMS Y [mm]",
            "RMS Plot"
        )
        
        beta_plot = create_plot(
            twi1[:,1],
            twi2[:,1],
            twi3[:,1],
            twi1[:,4],
            twi2[:,4],
            twi3[:,4],
            L"\beta_x",
            L"\beta_y",
            "Beta Functions"
        )
        
        alpha_plot = create_plot(
            twi1[:,2],
            twi2[:,2],
            twi3[:,2],
            twi1[:,5],
            twi2[:,5],
            twi3[:,5],
            L"\alpha_x",
            L"\alpha_y",
            "Alpha Functions"
        )
        
        emittance_plot = create_plot(
            twi1[:,3],
            twi2[:,3],
            twi3[:,3],
            twi1[:,6],
            twi2[:,6],
            twi3[:,6],
            L"\epsilon_x",
            L"\epsilon_y",
            "Emittance"
        )
        
        # Return all plots
        return Dict(
            "rms" => rms_plot,
            "beta" => beta_plot,
            "alpha" => alpha_plot,
            "emittance" => emittance_plot
        )
    else
        # Create one large figure
        f = Figure(size = (1200, 1200))
        
        # Create the main layout
        gl = f[1, 1] = GridLayout()
        legend_layout = f[1, 2] = GridLayout()
        
        # Define beam colors and labels
        beam_colors = [:blue, :orange, :green]
        beam_labels = [L"^{124}Xe^{49+}", L"^{124}Xe^{50+}", L"^{124}Xe^{51+}"]
        
        # Define element colors and heights
        element_colors = Dict(
            "SBEND" => :purple,
            "KQUAD" => :red,
            "ORBTRIM" => :green,
            "MARKER" => :gray,
            "DRIFT" => :white,
            "default" => :gray
        )
        
        element_heights = Dict(
            "SBEND" => 0.8,
            "KQUAD" => 0.7,
            "ORBTRIM" => 0.6,
            "MARKER" => 0.4,
            "DRIFT" => 0.3,
            "default" => 0.5
        )
        
        # Create plots in sequence (8 total + 1 lattice in the middle)
        # 1. RMS X
        ax_rms_x = Axis(gl[1, 1], xlabel = "", ylabel = "RMS X [mm]", title = "RMS Values")
        lines!(ax_rms_x, floor_length, beam1_rms[:,1] .* 1e3, color = beam_colors[1], linewidth = 2)
        lines!(ax_rms_x, floor_length, beam2_rms[:,1] .* 1e3, color = beam_colors[2], linewidth = 2)
        lines!(ax_rms_x, floor_length, beam3_rms[:,1] .* 1e3, color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_rms_x)
        
        # 2. RMS Y
        ax_rms_y = Axis(gl[2, 1], xlabel = "", ylabel = "RMS Y [mm]")
        lines!(ax_rms_y, floor_length, beam1_rms[:,3] .* 1e3, color = beam_colors[1], linewidth = 2)
        lines!(ax_rms_y, floor_length, beam2_rms[:,3] .* 1e3, color = beam_colors[2], linewidth = 2)
        lines!(ax_rms_y, floor_length, beam3_rms[:,3] .* 1e3, color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_rms_y)
        
        # 3. Beta X
        ax_beta_x = Axis(gl[3, 1], xlabel = "", ylabel = L"\beta_x", title = "Beta Functions")
        lines!(ax_beta_x, floor_length, twi1[:,1], color = beam_colors[1], linewidth = 2)
        lines!(ax_beta_x, floor_length, twi2[:,1], color = beam_colors[2], linewidth = 2)
        lines!(ax_beta_x, floor_length, twi3[:,1], color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_beta_x)
        
        # 4. Beta Y
        ax_beta_y = Axis(gl[4, 1], xlabel = "", ylabel = L"\beta_y")
        lines!(ax_beta_y, floor_length, twi1[:,4], color = beam_colors[1], linewidth = 2)
        lines!(ax_beta_y, floor_length, twi2[:,4], color = beam_colors[2], linewidth = 2)
        lines!(ax_beta_y, floor_length, twi3[:,4], color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_beta_y)
        
        # 5. Lattice plot in the middle
        floor_ax = Axis(gl[5, 1], 
                    xlabel = "",
                    ylabel = "",
                    yticklabelsvisible = false,
                    yticksvisible = false)
        
        # Track element types for legend
        element_types = Dict{String, Any}()
        
        # Draw the floorline
        curr_pos = 0.0
        for ele in lattice
            # Get element type
            ele_type = string(typeof(ele).name.name)
            
            # Determine color and height
            color = get(element_colors, ele_type, element_colors["default"])
            height = get(element_heights, ele_type, element_heights["default"])
            
            # Use a minimum width to prevent tiny elements and avoid overlap
            ele_width = ele.len #max(ele.len, 0.05)
            
            # Draw rectangle for element based on element type and properties
            if ele_type == "KQUAD" && hasfield(typeof(ele), :k1)
                # Use the original quad height, not drift height
                quad_height = element_heights["KQUAD"]
                if ele.k1 > 0
                    # Position touching x-axis and extending upward with full quad height
                    rect = Rect(curr_pos, 0, ele_width, quad_height)
                else
                    # Position touching x-axis and extending downward with full quad height
                    rect = Rect(curr_pos, -quad_height, ele_width, quad_height)
                end
            elseif ele_type == "KSEXT" && hasfield(typeof(ele), :k2)
                # Use the original sextupole height
                sext_height = get(element_heights, "KSEXT", element_heights["default"])
                if ele.k2 > 0
                    # Position touching x-axis and extending upward with full sextupole height
                    rect = Rect(curr_pos, 0, ele_width, sext_height)
                else
                    # Position touching x-axis and extending downward with full sextupole height
                    rect = Rect(curr_pos, -sext_height, ele_width, sext_height)
                end
            else
                # All other elements centered on x-axis as before
                rect = Rect(curr_pos, -height/2, ele_width, height)
            end
            
            element = poly!(floor_ax, rect, color = color, strokewidth = 1, strokecolor = :black)
            
            # Track for legend (only if not already tracked)
            if !haskey(element_types, ele_type)
                element_types[ele_type] = element
            end
            
            # Update position - use the element width to ensure no overlap
            curr_pos += ele_width
        end
        hidexdecorations!(floor_ax)
        
        # 6. Alpha X
        ax_alpha_x = Axis(gl[6, 1], xlabel = "", ylabel = L"\alpha_x", title = "Alpha Functions")
        lines!(ax_alpha_x, floor_length, twi1[:,2], color = beam_colors[1], linewidth = 2)
        lines!(ax_alpha_x, floor_length, twi2[:,2], color = beam_colors[2], linewidth = 2)
        lines!(ax_alpha_x, floor_length, twi3[:,2], color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_alpha_x)
        
        # 7. Alpha Y
        ax_alpha_y = Axis(gl[7, 1], xlabel = "", ylabel = L"\alpha_y")
        lines!(ax_alpha_y, floor_length, twi1[:,5], color = beam_colors[1], linewidth = 2)
        lines!(ax_alpha_y, floor_length, twi2[:,5], color = beam_colors[2], linewidth = 2)
        lines!(ax_alpha_y, floor_length, twi3[:,5], color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_alpha_y)
        
        # 8. Emittance X
        ax_emit_x = Axis(gl[8, 1], xlabel = "", ylabel = L"\epsilon_x", title = "Emittance")
        lines!(ax_emit_x, floor_length, twi1[:,3], color = beam_colors[1], linewidth = 2)
        lines!(ax_emit_x, floor_length, twi2[:,3], color = beam_colors[2], linewidth = 2)
        lines!(ax_emit_x, floor_length, twi3[:,3], color = beam_colors[3], linewidth = 2)
        hidexdecorations!(ax_emit_x)
        
        # 9. Emittance Y (bottom plot)
        ax_emit_y = Axis(gl[9, 1], xlabel = "z [m]", ylabel = L"\epsilon_y")
        lines!(ax_emit_y, floor_length, twi1[:,6], color = beam_colors[1], linewidth = 2)
        lines!(ax_emit_y, floor_length, twi2[:,6], color = beam_colors[2], linewidth = 2)
        lines!(ax_emit_y, floor_length, twi3[:,6], color = beam_colors[3], linewidth = 2)
        
        # Link all x axes
        for ax in [ax_rms_x, ax_rms_y, ax_beta_x, ax_beta_y, floor_ax, 
                ax_alpha_x, ax_alpha_y, ax_emit_x, ax_emit_y]
            linkxaxes!(ax, ax_rms_x)
        end
        
        # Add beam legend
        legend_entries = [
            LineElement(color = c, linewidth = 2) for c in beam_colors
        ]
        
        leg = Legend(legend_layout[1, 1], legend_entries, beam_labels, "Beam Types")
        
        # Add element legend (excluding MARKER and DRIFT)
        element_entries = []
        element_names = []
        sorted_names = sort(collect(keys(element_types)))
        for name in sorted_names
            if !(name in ["MARKER", "DRIFT"])
                element = element_types[name]
                push!(element_entries, element)
                push!(element_names, name)
            end
        end
        
        if !isempty(element_entries)
            element_leg = Legend(legend_layout[2, 1], element_entries, element_names, "Element Types")
        end
        
        # Adjust spacing between legends
        rowgap!(legend_layout, 5)
        
        # Set row sizes (make lattice row smaller)
        rowsize!(gl, 5, 30)  # Lattice plot is smaller
        
        
        return f
    end
end

# Create the lattice and beams
lattice, beam1, beam2, beam3 = create_lattice();

begin
    p1 = visualize_beam_properties(beam1);
    p2 = visualize_beam_properties(beam2);
    p3 = visualize_beam_properties(beam3);
    display(p1)
    display(p2)
    display(p3)
end;


# Run simulation
lattice, beam1, beam2, beam3, beam1_rms, beam2_rms, beam3_rms, 
twi1, twi2, twi3, floor_distance, expanded_indices = run_simulation();

# Plot using original indexing
end_ele = 200
plots = plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, false);
plots["rms"] 
plot_multibeam(end_ele, beam1_rms, beam2_rms, beam3_rms, twi1, twi2, twi3, floor_distance,lattice, expanded_indices, true)
