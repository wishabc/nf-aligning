// Workaround, so when we groupTuple later, 
// it knows how many objects in the group are going to be
def set_key_for_group_tuple(ch) {
  channel_1.groupTuple()
    .map(key, files -> tuple(groupKey(key, files.size()), files))
    .transpose()
}