
def get_container(file_name) {
  parent = file(file_name).parent
  old_parent = file(file_name).toRealPath().parent
  container = "--bind ${parent},${old_parent}"
}

// Workaround, so when we groupTuple later, 
// New key contains info on how many objects are in the group
def set_key_for_group_tuple(ch) {
  ch.groupTuple()
  .map{ it -> tuple(groupKey(it[0], it[1].size()), *it[1..(it.size()-1)]) }
  .transpose()
}

fastaContainer = get_container(params.genome_fasta_file)