# This file was created automatically by SWIG.
# Don't modify this file, modify the SWIG interface instead.
# This file is compatible with both classic and new-style classes.

import _pbs

def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "this"):
        if isinstance(value, class_type):
            self.__dict__[name] = value.this
            if hasattr(value,"thisown"): self.__dict__["thisown"] = value.thisown
            del value.thisown
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static) or hasattr(self,name) or (name == "thisown"):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError,name

import types
try:
    _object = types.ObjectType
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0
del types



new_attrl = _pbs.new_attrl

new_attropl = _pbs.new_attropl

new_batch_status = _pbs.new_batch_status

get_error = _pbs.get_error
ATTR_a = _pbs.ATTR_a
ATTR_c = _pbs.ATTR_c
ATTR_e = _pbs.ATTR_e
ATTR_g = _pbs.ATTR_g
ATTR_h = _pbs.ATTR_h
ATTR_j = _pbs.ATTR_j
ATTR_k = _pbs.ATTR_k
ATTR_l = _pbs.ATTR_l
ATTR_m = _pbs.ATTR_m
ATTR_o = _pbs.ATTR_o
ATTR_p = _pbs.ATTR_p
ATTR_q = _pbs.ATTR_q
ATTR_r = _pbs.ATTR_r
ATTR_u = _pbs.ATTR_u
ATTR_v = _pbs.ATTR_v
ATTR_A = _pbs.ATTR_A
ATTR_M = _pbs.ATTR_M
ATTR_N = _pbs.ATTR_N
ATTR_S = _pbs.ATTR_S
ATTR_depend = _pbs.ATTR_depend
ATTR_inter = _pbs.ATTR_inter
ATTR_stagein = _pbs.ATTR_stagein
ATTR_stageout = _pbs.ATTR_stageout
ATTR_ctime = _pbs.ATTR_ctime
ATTR_exechost = _pbs.ATTR_exechost
ATTR_mtime = _pbs.ATTR_mtime
ATTR_qtime = _pbs.ATTR_qtime
ATTR_session = _pbs.ATTR_session
ATTR_euser = _pbs.ATTR_euser
ATTR_egroup = _pbs.ATTR_egroup
ATTR_hashname = _pbs.ATTR_hashname
ATTR_hopcount = _pbs.ATTR_hopcount
ATTR_security = _pbs.ATTR_security
ATTR_sched_hint = _pbs.ATTR_sched_hint
ATTR_substate = _pbs.ATTR_substate
ATTR_name = _pbs.ATTR_name
ATTR_owner = _pbs.ATTR_owner
ATTR_used = _pbs.ATTR_used
ATTR_state = _pbs.ATTR_state
ATTR_queue = _pbs.ATTR_queue
ATTR_server = _pbs.ATTR_server
ATTR_maxrun = _pbs.ATTR_maxrun
ATTR_total = _pbs.ATTR_total
ATTR_comment = _pbs.ATTR_comment
ATTR_cookie = _pbs.ATTR_cookie
ATTR_qrank = _pbs.ATTR_qrank
ATTR_altid = _pbs.ATTR_altid
ATTR_etime = _pbs.ATTR_etime
ATTR_aclgren = _pbs.ATTR_aclgren
ATTR_aclgroup = _pbs.ATTR_aclgroup
ATTR_aclhten = _pbs.ATTR_aclhten
ATTR_aclhost = _pbs.ATTR_aclhost
ATTR_acluren = _pbs.ATTR_acluren
ATTR_acluser = _pbs.ATTR_acluser
ATTR_altrouter = _pbs.ATTR_altrouter
ATTR_chkptmin = _pbs.ATTR_chkptmin
ATTR_enable = _pbs.ATTR_enable
ATTR_fromroute = _pbs.ATTR_fromroute
ATTR_killdelay = _pbs.ATTR_killdelay
ATTR_maxgrprun = _pbs.ATTR_maxgrprun
ATTR_maxque = _pbs.ATTR_maxque
ATTR_maxuserrun = _pbs.ATTR_maxuserrun
ATTR_qtype = _pbs.ATTR_qtype
ATTR_rescassn = _pbs.ATTR_rescassn
ATTR_rescdflt = _pbs.ATTR_rescdflt
ATTR_rescmax = _pbs.ATTR_rescmax
ATTR_rescmin = _pbs.ATTR_rescmin
ATTR_rndzretry = _pbs.ATTR_rndzretry
ATTR_routedest = _pbs.ATTR_routedest
ATTR_routeheld = _pbs.ATTR_routeheld
ATTR_routewait = _pbs.ATTR_routewait
ATTR_routeretry = _pbs.ATTR_routeretry
ATTR_routelife = _pbs.ATTR_routelife
ATTR_rsvexpdt = _pbs.ATTR_rsvexpdt
ATTR_rsvsync = _pbs.ATTR_rsvsync
ATTR_start = _pbs.ATTR_start
ATTR_count = _pbs.ATTR_count
ATTR_number = _pbs.ATTR_number
ATTR_reqprop = _pbs.ATTR_reqprop
ATTR_aclroot = _pbs.ATTR_aclroot
ATTR_managers = _pbs.ATTR_managers
ATTR_dfltque = _pbs.ATTR_dfltque
ATTR_defnode = _pbs.ATTR_defnode
ATTR_locsvrs = _pbs.ATTR_locsvrs
ATTR_logevents = _pbs.ATTR_logevents
ATTR_logfile = _pbs.ATTR_logfile
ATTR_mailfrom = _pbs.ATTR_mailfrom
ATTR_nodepack = _pbs.ATTR_nodepack
ATTR_operators = _pbs.ATTR_operators
ATTR_queryother = _pbs.ATTR_queryother
ATTR_resccost = _pbs.ATTR_resccost
ATTR_rescavail = _pbs.ATTR_rescavail
ATTR_schedit = _pbs.ATTR_schedit
ATTR_scheduling = _pbs.ATTR_scheduling
ATTR_status = _pbs.ATTR_status
ATTR_syscost = _pbs.ATTR_syscost
ATTR_pingrate = _pbs.ATTR_pingrate
ATTR_ndchkrate = _pbs.ATTR_ndchkrate
ATTR_procpack = _pbs.ATTR_procpack
ATTR_NODE_state = _pbs.ATTR_NODE_state
ATTR_NODE_np = _pbs.ATTR_NODE_np
ATTR_NODE_properties = _pbs.ATTR_NODE_properties
ATTR_NODE_ntype = _pbs.ATTR_NODE_ntype
ATTR_NODE_jobs = _pbs.ATTR_NODE_jobs
CHECKPOINT_UNSPECIFIED = _pbs.CHECKPOINT_UNSPECIFIED
NO_HOLD = _pbs.NO_HOLD
NO_JOIN = _pbs.NO_JOIN
NO_KEEP = _pbs.NO_KEEP
MAIL_AT_ABORT = _pbs.MAIL_AT_ABORT
DELDELAY = _pbs.DELDELAY
USER_HOLD = _pbs.USER_HOLD
OTHER_HOLD = _pbs.OTHER_HOLD
SYSTEM_HOLD = _pbs.SYSTEM_HOLD
ND_free = _pbs.ND_free
ND_offline = _pbs.ND_offline
ND_down = _pbs.ND_down
ND_reserve = _pbs.ND_reserve
ND_job_exclusive = _pbs.ND_job_exclusive
ND_job_sharing = _pbs.ND_job_sharing
ND_busy = _pbs.ND_busy
ND_state_unknown = _pbs.ND_state_unknown
ND_timeshared = _pbs.ND_timeshared
ND_cluster = _pbs.ND_cluster
MAX_ENCODE_BFR = _pbs.MAX_ENCODE_BFR
MGR_CMD_CREATE = _pbs.MGR_CMD_CREATE
MGR_CMD_DELETE = _pbs.MGR_CMD_DELETE
MGR_CMD_SET = _pbs.MGR_CMD_SET
MGR_CMD_UNSET = _pbs.MGR_CMD_UNSET
MGR_CMD_LIST = _pbs.MGR_CMD_LIST
MGR_CMD_PRINT = _pbs.MGR_CMD_PRINT
MGR_CMD_ACTIVE = _pbs.MGR_CMD_ACTIVE
MGR_OBJ_NONE = _pbs.MGR_OBJ_NONE
MGR_OBJ_SERVER = _pbs.MGR_OBJ_SERVER
MGR_OBJ_QUEUE = _pbs.MGR_OBJ_QUEUE
MGR_OBJ_JOB = _pbs.MGR_OBJ_JOB
MGR_OBJ_NODE = _pbs.MGR_OBJ_NODE
MSG_OUT = _pbs.MSG_OUT
MSG_ERR = _pbs.MSG_ERR
SHUT_SIG = _pbs.SHUT_SIG
SHUT_IMMEDIATE = _pbs.SHUT_IMMEDIATE
SHUT_DELAY = _pbs.SHUT_DELAY
SHUT_QUICK = _pbs.SHUT_QUICK
SIG_RESUME = _pbs.SIG_RESUME
SIG_SUSPEND = _pbs.SIG_SUSPEND
PBS_MAXHOSTNAME = _pbs.PBS_MAXHOSTNAME
MAXPATHLEN = _pbs.MAXPATHLEN
MAXNAMLEN = _pbs.MAXNAMLEN
PBS_MAXUSER = _pbs.PBS_MAXUSER
PBS_MAXGRPN = _pbs.PBS_MAXGRPN
PBS_MAXQUEUENAME = _pbs.PBS_MAXQUEUENAME
PBS_MAXSERVERNAME = _pbs.PBS_MAXSERVERNAME
PBS_MAXSEQNUM = _pbs.PBS_MAXSEQNUM
PBS_MAXPORTNUM = _pbs.PBS_MAXPORTNUM
PBS_MAXSVRJOBID = _pbs.PBS_MAXSVRJOBID
PBS_MAXCLTJOBID = _pbs.PBS_MAXCLTJOBID
PBS_MAXDEST = _pbs.PBS_MAXDEST
PBS_MAXROUTEDEST = _pbs.PBS_MAXROUTEDEST
PBS_USE_IFF = _pbs.PBS_USE_IFF
PBS_INTERACTIVE = _pbs.PBS_INTERACTIVE
PBS_TERM_BUF_SZ = _pbs.PBS_TERM_BUF_SZ
PBS_TERM_CCA = _pbs.PBS_TERM_CCA
PBS_BATCH_SERVICE_NAME = _pbs.PBS_BATCH_SERVICE_NAME
PBS_BATCH_SERVICE_PORT = _pbs.PBS_BATCH_SERVICE_PORT
PBS_BATCH_SERVICE_NAME_DIS = _pbs.PBS_BATCH_SERVICE_NAME_DIS
PBS_BATCH_SERVICE_PORT_DIS = _pbs.PBS_BATCH_SERVICE_PORT_DIS
PBS_MOM_SERVICE_NAME = _pbs.PBS_MOM_SERVICE_NAME
PBS_MOM_SERVICE_PORT = _pbs.PBS_MOM_SERVICE_PORT
PBS_MANAGER_SERVICE_NAME = _pbs.PBS_MANAGER_SERVICE_NAME
PBS_MANAGER_SERVICE_PORT = _pbs.PBS_MANAGER_SERVICE_PORT
PBS_SCHEDULER_SERVICE_NAME = _pbs.PBS_SCHEDULER_SERVICE_NAME
PBS_SCHEDULER_SERVICE_PORT = _pbs.PBS_SCHEDULER_SERVICE_PORT
SET = _pbs.SET
UNSET = _pbs.UNSET
INCR = _pbs.INCR
DECR = _pbs.DECR
EQ = _pbs.EQ
NE = _pbs.NE
GE = _pbs.GE
GT = _pbs.GT
LE = _pbs.LE
LT = _pbs.LT
DFLT = _pbs.DFLT
class attrl(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, attrl, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, attrl, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C attrl instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["name"] = _pbs.attrl_name_set
    __swig_getmethods__["name"] = _pbs.attrl_name_get
    if _newclass:name = property(_pbs.attrl_name_get, _pbs.attrl_name_set)
    __swig_setmethods__["resource"] = _pbs.attrl_resource_set
    __swig_getmethods__["resource"] = _pbs.attrl_resource_get
    if _newclass:resource = property(_pbs.attrl_resource_get, _pbs.attrl_resource_set)
    __swig_setmethods__["value"] = _pbs.attrl_value_set
    __swig_getmethods__["value"] = _pbs.attrl_value_get
    if _newclass:value = property(_pbs.attrl_value_get, _pbs.attrl_value_set)
    __swig_setmethods__["op"] = _pbs.attrl_op_set
    __swig_getmethods__["op"] = _pbs.attrl_op_get
    if _newclass:op = property(_pbs.attrl_op_get, _pbs.attrl_op_set)
    def __str__(*args): return _pbs.attrl___str__(*args)

class attrlPtr(attrl):
    def __init__(self, this):
        _swig_setattr(self, attrl, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, attrl, 'thisown', 0)
        _swig_setattr(self, attrl,self.__class__,attrl)
_pbs.attrl_swigregister(attrlPtr)

class attropl(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, attropl, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, attropl, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C attropl instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["name"] = _pbs.attropl_name_set
    __swig_getmethods__["name"] = _pbs.attropl_name_get
    if _newclass:name = property(_pbs.attropl_name_get, _pbs.attropl_name_set)
    __swig_setmethods__["resource"] = _pbs.attropl_resource_set
    __swig_getmethods__["resource"] = _pbs.attropl_resource_get
    if _newclass:resource = property(_pbs.attropl_resource_get, _pbs.attropl_resource_set)
    __swig_setmethods__["value"] = _pbs.attropl_value_set
    __swig_getmethods__["value"] = _pbs.attropl_value_get
    if _newclass:value = property(_pbs.attropl_value_get, _pbs.attropl_value_set)
    __swig_setmethods__["op"] = _pbs.attropl_op_set
    __swig_getmethods__["op"] = _pbs.attropl_op_get
    if _newclass:op = property(_pbs.attropl_op_get, _pbs.attropl_op_set)
    def __str__(*args): return _pbs.attropl___str__(*args)

class attroplPtr(attropl):
    def __init__(self, this):
        _swig_setattr(self, attropl, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, attropl, 'thisown', 0)
        _swig_setattr(self, attropl,self.__class__,attropl)
_pbs.attropl_swigregister(attroplPtr)

class batch_status(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, batch_status, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, batch_status, name)
    def __init__(self): raise RuntimeError, "No constructor defined"
    def __repr__(self):
        return "<%s.%s; proxy of C batch_status instance at %s>" % (self.__class__.__module__, self.__class__.__name__, self.this,)
    __swig_setmethods__["name"] = _pbs.batch_status_name_set
    __swig_getmethods__["name"] = _pbs.batch_status_name_get
    if _newclass:name = property(_pbs.batch_status_name_get, _pbs.batch_status_name_set)
    __swig_setmethods__["attribs"] = _pbs.batch_status_attribs_set
    __swig_getmethods__["attribs"] = _pbs.batch_status_attribs_get
    if _newclass:attribs = property(_pbs.batch_status_attribs_get, _pbs.batch_status_attribs_set)
    __swig_setmethods__["text"] = _pbs.batch_status_text_set
    __swig_getmethods__["text"] = _pbs.batch_status_text_get
    if _newclass:text = property(_pbs.batch_status_text_get, _pbs.batch_status_text_set)

class batch_statusPtr(batch_status):
    def __init__(self, this):
        _swig_setattr(self, batch_status, 'this', this)
        if not hasattr(self,"thisown"): _swig_setattr(self, batch_status, 'thisown', 0)
        _swig_setattr(self, batch_status,self.__class__,batch_status)
_pbs.batch_status_swigregister(batch_statusPtr)

RESOURCE_T_NULL = _pbs.RESOURCE_T_NULL
RESOURCE_T_ALL = _pbs.RESOURCE_T_ALL

avail = _pbs.avail

pbs_asyrunjob = _pbs.pbs_asyrunjob

pbs_alterjob = _pbs.pbs_alterjob

pbs_connect = _pbs.pbs_connect

pbs_query_max_connections = _pbs.pbs_query_max_connections

pbs_default = _pbs.pbs_default

pbs_deljob = _pbs.pbs_deljob

pbs_disconnect = _pbs.pbs_disconnect

pbs_holdjob = _pbs.pbs_holdjob

pbs_locjob = _pbs.pbs_locjob

pbs_manager = _pbs.pbs_manager

pbs_movejob = _pbs.pbs_movejob

pbs_msgjob = _pbs.pbs_msgjob

pbs_orderjob = _pbs.pbs_orderjob

pbs_rescquery = _pbs.pbs_rescquery

pbs_rescreserve = _pbs.pbs_rescreserve

pbs_rescrelease = _pbs.pbs_rescrelease

pbs_rerunjob = _pbs.pbs_rerunjob

pbs_rlsjob = _pbs.pbs_rlsjob

pbs_runjob = _pbs.pbs_runjob

pbs_selectjob = _pbs.pbs_selectjob

pbs_sigjob = _pbs.pbs_sigjob

pbs_statfree = _pbs.pbs_statfree

pbs_statjob = _pbs.pbs_statjob

pbs_selstat = _pbs.pbs_selstat

pbs_statque = _pbs.pbs_statque

pbs_statserver = _pbs.pbs_statserver

pbs_statnode = _pbs.pbs_statnode

pbs_submit = _pbs.pbs_submit

pbs_terminate = _pbs.pbs_terminate

totpool = _pbs.totpool

usepool = _pbs.usepool

openrm = _pbs.openrm

closerm = _pbs.closerm

downrm = _pbs.downrm

configrm = _pbs.configrm

addreq = _pbs.addreq

allreq = _pbs.allreq

flushreq = _pbs.flushreq

activereq = _pbs.activereq

fullresp = _pbs.fullresp

getreq = _pbs.getreq
LOG_BUF_SIZE = _pbs.LOG_BUF_SIZE

log_close = _pbs.log_close

log_err = _pbs.log_err

log_event = _pbs.log_event

log_open = _pbs.log_open

log_record = _pbs.log_record

setup_env = _pbs.setup_env

chk_file_sec = _pbs.chk_file_sec
PBSEVENT_ERROR = _pbs.PBSEVENT_ERROR
PBSEVENT_SYSTEM = _pbs.PBSEVENT_SYSTEM
PBSEVENT_ADMIN = _pbs.PBSEVENT_ADMIN
PBSEVENT_JOB = _pbs.PBSEVENT_JOB
PBSEVENT_JOB_USAGE = _pbs.PBSEVENT_JOB_USAGE
PBSEVENT_SECURITY = _pbs.PBSEVENT_SECURITY
PBSEVENT_SCHED = _pbs.PBSEVENT_SCHED
PBSEVENT_DEBUG = _pbs.PBSEVENT_DEBUG
PBSEVENT_DEBUG2 = _pbs.PBSEVENT_DEBUG2
PBSEVENT_FORCE = _pbs.PBSEVENT_FORCE
PBS_EVENTCLASS_SERVER = _pbs.PBS_EVENTCLASS_SERVER
PBS_EVENTCLASS_QUEUE = _pbs.PBS_EVENTCLASS_QUEUE
PBS_EVENTCLASS_JOB = _pbs.PBS_EVENTCLASS_JOB
PBS_EVENTCLASS_REQUEST = _pbs.PBS_EVENTCLASS_REQUEST
PBS_EVENTCLASS_FILE = _pbs.PBS_EVENTCLASS_FILE
PBS_EVENTCLASS_ACCT = _pbs.PBS_EVENTCLASS_ACCT
PBS_EVENTCLASS_NODE = _pbs.PBS_EVENTCLASS_NODE
PBSEVENT_MASK = _pbs.PBSEVENT_MASK
#  PBS python interface
#  Author: Bas van der Vlies <basv@sara.nl>
#  Date  : 27 Feb 2002
#  Desc. : This is python wrapper class for getting the resource
#          mom values.
#
# CVS info
# $Id: pbs.py 144 2006-11-16 13:35:53Z bas $
# $Date: 2002/10/21 14:14:47 $
# $Revision: 1.6 $
#
import string
import types

# Default linux resources to get from the mom
#
default_linux_res = [   
	"availmem",	# available memory size in KB
	"ideal_load",	# static ideal_load value
	"loadave",      # the current load average
	"max_load",	# static max_load value
	"ncpus",        # number of cpus
	"physmem",      # physical memory size in KB
	"resi",		# resident memory size for a pid or session in KB
	"totmem",	# total memory size in KB
	"walltime",	# wall clock time for a pid
]

# Default irix6 resources to get from the mom
#
default_irix6_res = [   
	"availmem",	# available memory size in KB
	"loadave",      # the current load average
	"ncpus",        # number of cpus
	"physmem",      # physical memory size in KB
	"resi",		# resident memory size for a pid or session in KB
	"walltime",	# wall clock time for a pid
	"quota",	# quota information (sizes in KB)
]

default_mom_res = [   
	"arch",		# the architecture of the machine
	"uname",	# the architecture of the machine
        "cput",		# cpu time for a pid or session
	"idletime",	# seconds of idle time
	"mem",		# memory size for a pid or session in KB
	"sessions",	# list of sessions in the system
	"pids",         # list of pids in a session
	"nsessions",	# number of sessions in the system
	"nusers",	# number of users in the system
	"size",		# size of a file or filesystem
	"host",		# Name  of host on which job should be run 
	"nodes",	# Number and/or type of nodes to be reserved for exclusive use by the job
	"other",	# Allows a  user  to  specify  site  specific  information
	"software",	# Allows a user to specify software required by the job
]

def check_resp(dict, str):
  """
  Check the daemon response. If we have no permission to
  query the values then we got a 'None' response. Else
  if we supplied a keyword that does not exits we get a
  '?' response
  """
  if not str:
    return
      
  key, val = string.split(str, '=')
  key = string.strip(key)
  val = string.strip(val)

  # Did we got a valid response
  #
  if not val[0] == '?':
    dict[key] = val

def use_default_keywords(id, d):
  """
  Get the default values from the mom daemon
  """
  for res in default_mom_res:
    addreq(id, res)
    resp = getreq(id)
    check_resp(d, resp)

  # Do not proceed if we have an empty dictionary
  #
  if not d:
    return

  if d['arch' ] == 'linux':
    for res in default_linux_res:
      addreq(id, res)
      resp = getreq(id)
      check_resp(d, resp)

def use_user_keywords(id, d, l):
  for res in l:
    if type(res) is types.StringType:
      addreq(id, res)
      resp = getreq(id)
      check_resp(d, resp)
    else:
      raise TypeError, 'Expected a string got %s :%s' %(type(res), res) 

def get_mom_values(id, list = None):
  """
  This function will query the mom with a default resmon keywords
  and 'arch' depended keywords. Supported archs are:
    linux
    irix6
  User can also supply their own list of keywords as second parameter.
  arguments:
    id   : connection number with mom daemon on a node
    list : optional parameter. If supplied then use this. A list
           of mom keywords.
  """

  d = {}
  if not list:
    use_default_keywords(id, d)
  else:
    use_user_keywords(id, d , list)
     
  return d

def version():
  """
  Returns the pbs python interface version as a string. 
  """
  return '2.9.4'

# A useful dict with error codes to text
#
# SVN Info:
#	$Id: pbs.py 144 2006-11-16 13:35:53Z bas $
#
errors_txt = { 
	0 : 'no error',
	15001 :	 'Unknown Job Identifier',
	15002 : 'Undefined Attribute',
	15003 : 'attempt to set READ ONLY attribute',
	15004 : 'Invalid request',
	15005 : 'Unknown batch request',
	15006 : 'Too many submit retries',
	15007 : 'No permission',
	15008 : 'access from host not allowed',
	15009 : 'job already exists',
	15010 : 'system error occurred',
	15011 : 'internal server error occurred',
	15012 : 'parent job of dependent in rte que',
	15013 : 'unknown signal name',
	15014 : 'bad attribute value',
	15015 : 'Cannot modify attrib in run state',
	15016 : 'request invalid for job state',
	15018 : 'Unknown queue name',
	15019 : 'Invalid Credential in request',
	15020 : 'Expired Credential in request',
	15021 : 'Queue not enabled',
	15022 : 'No access permission for queue',
	15023 : 'Bad user - no password entry',
	15024 : 'Max hop count exceeded',
	15025 : 'Queue already exists',
	15026 : 'incompatable queue attribute type',
	15027 : 'Queue Busy (not empty)',
	15028 : 'Queue name too long',
	15029 : 'Feature',
	15030 : 'Cannot enable queue,needs add def',
	15031 : 'Protocol (ASN.1) error',
	15032 : 'Bad attribute list structure',
	15033 : 'No free connections',
	15034 : 'No server to connect to',
	15035 : 'Unknown resource',
	15036 : 'Job exceeds Queue resource limits',
	15037 : 'No Default Queue Defined',
	15038 : 'Job Not Rerunnable',
	15039 : 'Route rejected by all destinations',
	15040 : 'Time in Route Queue Expired',
	15041 : 'Request to MOM failed',
	15042 : '(qsub) cannot access script file',
	15043 : 'Stage In of files failed',
	15044 : 'Resources temporarily unavailable',
	15045 : 'Bad Group specified',
	15046 : 'Max number of jobs in queue',
	15047 : 'Checkpoint Busy, may be retries',
	15048 : 'Limit exceeds allowable',
	15049 : 'Bad Account attribute value',
	15050 : 'Job already in exit state',
	15051 : 'Job files not copied',
	15052 : 'unknown job id after clean init',
	15053 : 'No Master in Sync Set',
	15054 : 'Invalid dependency',
	15055 : 'Duplicate entry in List',
	15056 : 'Bad DIS based Request Protocol',
	15057 : 'cannot execute there',
	15058 : 'sister rejected',
	15059 : 'sister could not communicate',
	15060 : 'req rejected -server shutting down',
	15061 : 'not all tasks could checkpoint',
	15062 : 'Named node is not in the list',
	15063 : 'node-attribute not recognized',
	15064 : 'Server has no node list',
	15065 : 'Node name is too big',
	15066 : 'Node name already exists',
	15067 : 'Bad node-attribute value',
	15068 : 'State values are mutually exclusive',
	15069 : 'Error(s) during global modification of nodes',
	15070 : 'could not contact Mom',
	15071 : 'no time-shared nodes',
	15201 : 'resource unknown',
	15202 : 'parameter could not be used',
	15203 : 'a parameter needed did not exist',
	15204 : "something specified didn't exist",
	15205 : 'a system error occured',
	15206 : 'only part of reservation made'
}

def error():
  """
  Check if there is an error, if so fetch the error message string. 
  It says more then a number!
  """
  e = get_error()
  if errors_txt.has_key(e):
     return (e, errors_txt[e])
  else:
     return (e, "Could not find a text for this error, uhhh")

cvar = _pbs.cvar

