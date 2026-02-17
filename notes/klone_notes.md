# klone notes

## These are notes on the details of our klone HPC computing resources

These are generally available only to members of the UW Coastal Modeling Group. They are administered by Parker MacCready and eventually by Kate Hewett.

### Background
- klone HPC compute resources are sold in "slices", each of which has 32 cores.
- a physical node had 6 slices or 192 cores. We often request "exclusive" use of a node and use all the cores for a typical model run in the LiveOcean system.
- Each slice comes with 1 TB of disk space on klone.
- The cost of a slice is currently (2/2026) $4,288.44. You also generally have to pay "slot fees" of $1k per **node** per year (about $6k per year for our group), as per agreement with the College of the Environment, which covers part of the slot fees.
- Both the up-front costs and the slot fees have 0% IDC.
- We generally own "cpu-g2" slices (AMD), and they are good for 5 years. After that you have to buy them again.
- IMPORTANT: You cannot just buy more slices when you want; your unit (CoEnv in this case) has to be part of the hyak system and is only eligible for a set number of slices. Our ability to have as many slices as we do is because (i) we are grandfathered into the system because of purchases beginning with the start of hyak, (ii) because we paid for some of the full slot fees ourselves, and (iii) because of a separate agreement between CoEnv and Oceanography - which resulted in the "coenv" group below.
- So if we want to maintain our slice count we need to replace them when they expire.

---

### Table of Current Resources

| Group | cpu type | Slices (cpu) | Start-End Dates | Slot Fees & Notes |
| --- | --- | --- | --- | --- |
| coenv | cpu-g2| 10 (320) | 1/2025-1/21/2030 | $2k |
| coenv | cpu-g2| 2 (64) | 1/2026-1/2031 | (included in above) |
| macc | cpu-g2| 10 (320) | 1/2025-1/21/2030 | $0k (I paid slot fees up-front) |
| macc | cpu-g2| 5 (160) | 6/2024-6/13/2029 | $1k |
| macc | cpu-g2| 10 (320) | 1/2026-1/2031 | $2k |
| macc | compute | 5 (200) | 5/2023-5/15/2027 | $1k (REPLACE THESE!) |

So the total annual slot fees are $6k, I estimate.

---

### Thinking of these resources from a computing standpoint
- You have to have a job reside in a "Group".
- Assume a job takes one full node = 6 slices = 192 cpu
- in the coenv group we have 12 slices = 2 nodes
- in the macc group we have 25 slices = 4 nodes + one extra slice (all cpu-g2)
- we also have 5 of the old compute nodes, but these should be replaced as soon as we can, like May 2027.
- All told, we have resource for 6 full-node jobs to be running at the same time!
- In terms of disk space on klone I think we have 30 TB in the macc group, and 12 TB (?) in the coenv group.
- There is much more disk space (300 TB) for our group in kopah. See LO/Notes/kopah_notes.md for more info.
