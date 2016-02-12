from django.shortcuts import render
from django.http import HttpResponse

import json
import uuid
import logging
import StringIO

from pyVEP import VEP

from MutationInfo import MutationInfo 

print 'Initializing MutationInfo..'
mi = MutationInfo()
print 'Done'

root_logger = logging.getLogger()
root_logger.setLevel(logging.DEBUG)


# Create your views here.

def do_MutationInfo(request):
    '''
    test_ret = {'chrom': '12', 'source': 'NC_transcript', 'genome': 'GRCh37.p13', 'offset': 21355487, 'alt': 'G', 'ref': 'T'}
    '''

    print request.method

    if request.method == u'GET':

        return render(request, 'MutationInfoApp/MutationInfoApp.html', {

    	})

    elif request.method == u'POST':
        POST = request.POST

        variant_text = POST.getlist(u'variantText')[0]

#        variant_text = u'NM_006446.4:c.1198T>G'

        #Set a new logger
        string_buffer = StringIO.StringIO()
        logging_handler = logging.StreamHandler(string_buffer)
        logging_handler.setLevel(logging.DEBUG)
        this_uuid = str(uuid.uuid4())
        #formatter = logging.Formatter(this_uuid + ' - %(asctime)s - %(name)s - %(levelname)s - %(message)s')
        formatter = logging.Formatter(this_uuid + '%(levelname)s - %(message)s')
        logging_handler.setFormatter(formatter)
        root_logger.addHandler(logging_handler)

        #Run Mutation Info
        exception_message = 'MutationInfo run out of options!'
        mi_ret = None
        try:
            mi_ret = mi.get_info(variant_text)
        except Exception as e:
            exception_message = 'MutationInfo failed: ' + str(e)

        # Run VEP
        VEP_transcript_consequences = 'Transcript consequences not found'
        if type(mi_ret) is dict:
            # MutationInfo seems to run fine.
            # Try VEP
            VEP_variant = '%s:g.%s%s>%s' % tuple(map(str, (mi_ret['chrom'], mi_ret['offset'], mi_ret['ref'], mi_ret['alt'], )))
            #print 'VEP_variant:', VEP_variant
            VEP_ret = VEP(VEP_variant, mi_ret['genome'])
            VEP_msc = 'Could not run Variant Effect Predictor'
            if type(VEP_ret) is list:
                if len(VEP_ret) > 0:
                    if type(VEP_ret[0]) is dict:
                        if u'most_severe_consequence' in VEP_ret[0]:
                            VEP_msc = VEP_ret[0][u'most_severe_consequence']
                        if u'transcript_consequences' in VEP_ret[0]:
                            print VEP_ret[0]
                            VEP_transcript_consequences = [x for x in VEP_ret[0][u'transcript_consequences'] if (u'polyphen_score' in x) or (u'sift_score' in x)]
            elif type(VEP_ret) is dict:
                if u'error' in VEP_ret:
                    VEP_msc = VEP_ret[u'error']
            #END OF VEP


        #Get log messages
        log_messages = string_buffer.getvalue()
        log_messages = '\n'.join([x.replace(this_uuid, '') for x in log_messages.split('\n') if this_uuid in x])
        root_logger.removeHandler(logging_handler)

        if type(mi_ret) is dict:
            mi_ret['VEP_msc'] = VEP_msc
            #print VEP_transcript_consequences
            mi_ret['VEP_transcript_consequences'] = json.dumps(VEP_transcript_consequences, indent=4)
            #print mi_ret['VEP_transcript_consequences']
            mi_ret['log'] = log_messages
            mi_ret['success'] = True

        elif mi_ret is None:
            mi_ret = {
                'log': log_messages,
                'error_msg': exception_message,
                'success': False,
            }

        to_ret_json = json.dumps(mi_ret)
        return HttpResponse(to_ret_json, content_type='application/json')



