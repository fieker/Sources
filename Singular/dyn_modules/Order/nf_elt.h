#ifndef NF_ELT_HPP
#define NF_ELT_HPP

extern n_coeffType nf_type;

typedef struct {
  number den; //of type n_Z
  number elt; //of type n_NFord
} nf_elt_t;

#define nf_elt_den(a) ((nf_elt_t*)a->den)
#define nf_elt_num(a) ((nf_elt_t*)a->num)

BOOLEAN n_nfInit(coeffs r,  void * parameter);

#endif
