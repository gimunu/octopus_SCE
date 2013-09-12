#ifdef TEMPLATE_H
#undef ITEMPLATE_NAME
#undef TEMPLATE
#else
#define TEMPLATE_H
#define CONCAT(T,N) T ## _ ## N
#define EVAL(T,N) CONCAT(T,N)
#endif

#ifdef TEMPLATE_NAME
#ifdef TEMPLATE_TYPE
#define ITEMPLATE_NAME EVAL(TEMPLATE_NAME,TEMPLATE_TYPE)
#else
#define ITEMPLATE_NAME TEMPLATE_NAME
#endif
#else
#ifdef TEMPLATE_TYPE
#define ITEMPLATE_NAME TEMPLATE_TYPE
#else
#define ITEMPLATE_NAME
#endif
#endif
#define TEMPLATE(NAME) EVAL(ITEMPLATE_NAME,NAME)

#ifdef SUBTEMPLATE_NAME
#ifdef SUBTEMPLATE_TYPE
#define ISUBTEMPLATE_NAME EVAL(SUBTEMPLATE_NAME,SUBTEMPLATE_TYPE)
#else
#define ISUBTEMPLATE_NAME SUBTEMPLATE_NAME
#endif
#else
#define ISUBTEMPLATE_NAME null
#endif
#define SUBTEMPLATE(NAME) EVAL(ISUBTEMPLATE_NAME,NAME)

!! Local Variables:
!! mode: f90
!! End:

