# LO Design Notes

The LiveOcean modeling framework is designed to meet a number of goals:
- Run daily forecasts
- Run long hindcasts
- Run idealized cases
- Run pre- an post-processing jobs

For the Coastal Modeling Group at the University of Washington these jobs are usually associated with simulation of the coastal and estuarine waters near Washington State, but they can be applied anywhere. They are designed around the ROMS (Regional Ocean Modeling System) community ocean model.

To meet our goals a lot of thought has gone into code design. Here I will try to describe this design in the hope that it will guide future scientists who are embarking on projects of similar complexity.

---

The primary software architecture pattern is "pipe and filter," like linux. Each important task is conceptualized as:
```
Input => Code => Output
```
- Input may come from external web sources, our data files (kept in LO_data), ROMS output (kept in LO_roms), or the output of previous tasks (kept in LO_output).
- Code (kept in LO) is the only part of the system stored in GitHub.
- Output is typically binary files that go to LO_output.

To do any big jobs you have to string together many of these "important tasks". We keep tight control on the format of input and output so that changes to the code do not break the system.

---
