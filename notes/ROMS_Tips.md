# ROMS Tips

## This is a place for miscellaneous useful information about ROMS

#### Code is in an svn repo called LO_roms_source, in the folder ROMS, unless otherwise noted

---

`Utility/dateclock.F` gives info about acceptable naming of time units. In particular see the subroutine "time_units" at line 1345, and its comment section for "valid unit attributes".

`External/varinfo.yaml` gives the required names, dimensions, units, etc. for all input and output variables.

`LO_roms_source_alt/varinfo/varinfo.yaml` is the slightly edited version of the above that we use in the LO system.  It is used both for creating forcing and for running ROMS.

`Include/cppdefs.h` gives the definitions of all the possible cpp flags to define in the LO_roms_user/[ex_name]/[ex_name].h file you use for compiling.

`Include/globaldefs.h` is where defaults and overrides of the cpp flags are handled. This is a great place to look if some flag is not doing what you think it should.

Question: how do you know what order ROMS will try reading forcing files?

Question: how do you know what forcing variables ROMS will look for?

Answer: look in varinfo.yaml for the "index_code" of a variable you are interested in.  E.g. for river transport it is idRtra. Then use grep to see where this is referenced in the code:
```
(loenv) MacBook-Pro-2:ROMS pm8$ grep -rnI idRtra *
Adjoint/ad_get_data.F:107:        CALL get_ngfldr (ng, iADM, idRtra, SSF(ng)%ncid,                &
Adjoint/ad_set_data.F:708:          CALL set_ngfldr (ng, iADM, idRtra, 1, Nsrc(ng), 1,            &
External/varinfo.yaml:1843:    index_code:     idRtra
Modules/mod_ncparam.F:261:      integer  :: idRtra        ! river runoff mass transport
Modules/mod_ncparam.F:1668:          CASE ('idRtra')
Modules/mod_ncparam.F:1669:            idRtra=varid
Nonlinear/get_data.F:102:        CALL get_ngfld (ng, iNLM, idRtra, SSF(ng)%ncid,                 &
Nonlinear/set_data.F:130:          CALL set_ngfld (ng, iNLM, idRtra, 1, Nsrc(ng), 1,             &
Representer/rp_get_data.F:104:        CALL get_ngfld (ng, iRPM, idRtra, SSF(ng)%ncid,                 &
Representer/rp_set_data.F:812:          CALL set_ngfld (ng, iRPM, idRtra, 1, Nsrc(ng), 1,             &
Tangent/tl_get_data.F:104:        CALL get_ngfld (ng, iTLM, idRtra, SSF(ng)%ncid,                 &
Tangent/tl_set_data.F:713:          CALL set_ngfld (ng, iTLM, idRtra, 1, Nsrc(ng), 1,             &
```
The likely one to look at in this list is `Nonlinear/get_data.F`, and this will give the answers about exactly what variables are expected and in what order (e.g. point sinks/sources and then forcing). An even more concise version of this would be on klone in `LO_roms_user/[ex_name]/Build_roms/get_data.f90`, which is get_data.F after it has been parsed using the cpp flags.
