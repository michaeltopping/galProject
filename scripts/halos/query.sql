select DES.haloId as descendand_id,
    DES.np as descendant_mass,
    PROG.*
from millimil..MPAHalo DES,
    millimil..MPAHalo PROG
where DES.snapnum = 63
  and DES.np > 1e5
  and PROG.haloId between DES.haloId and DES.lastprogenitorId
  and PROG.redshift between 3.05 and 3.12
  and PROG.np > 100
order by DES.np desc, PROG.np desc
