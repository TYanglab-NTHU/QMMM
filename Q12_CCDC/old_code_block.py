
# add the neighborhood of the functional group atoms
# build a list of all atoms neighboring the functional groups

qm_cnt += 1
new_at = []
for fgat in shell03:
  for at in fgat.neighbours:
    if at not in shell02:
      new_at.append(at)
      if len(at.rings) != 0:
        if at.rings[0].is_aromatic:
          for z in at.rings[0].atoms:
            new_at.append(z)
if len(new_at) > 0:
  new_at=list(set(new_at))
  shell03 += new_at
  shell03=list(set(shell03))
  del new_at
fname = "qm_at%02i.xyz" % qm_cnt
list_xyz(shell03, fname, "Center, func. groups, neighbors")
# merge the results
qm += shell03
qm=list(set(qm))
