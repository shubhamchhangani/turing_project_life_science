"""Read short, noncredential progress summaries from the isolated task trials."""
import subprocess
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("--prefix",default="inverse-receptor-cooperativity__")
args=parser.parse_args()

reader=r'''
import json
from pathlib import Path
for filename in ['codex.txt','claude-code.txt']:
    p=Path('/logs/agent')/filename
    if not p.is_file():continue
    lines=p.read_text(errors='replace').splitlines()
    events=[]
    for line in lines:
        try:events.append(json.loads(line))
        except ValueError:pass
    print(filename,'events:',len(events))
    for event in events:
        if event.get('type')=='system' and event.get('subtype')=='init':
            print('Model:',event.get('model'))
    summaries=[]
    def show(*parts):summaries.append(' '.join(map(str,parts)))
    for event in events:
        typ=event.get('type')
        if filename=='codex.txt':
            item=event.get('item',{})
            if item.get('type')=='command_execution':
                if typ=='item.completed':show('command',item.get('command','')[:150], 'exit',item.get('exit_code'))
            elif item.get('text'):show(typ,item['text'][-500:])
            elif typ in ['error','turn.failed']:show(typ,event.get('error',event.get('message')))
        elif typ=='assistant':
            for item in event.get('message',{}).get('content',[]):
                if item.get('type')=='text':show('assistant',item.get('text','')[-500:])
                elif item.get('type')=='tool_use':
                    show('tool',item.get('name'),str(item.get('input',{}).get('command',item.get('input',{}).get('file_path','')))[:150])
        elif typ=='result':show('result',event.get('subtype'),str(event.get('result',''))[-500:])
    for summary in summaries[-4:]:print(summary)
'''
result=subprocess.run(['docker','ps','--format','{{.Names}}'],capture_output=True,text=True,check=True)
names=[n for n in result.stdout.splitlines() if n.startswith(args.prefix)]
for name in names:
    print(name,flush=True)
    run=subprocess.run(['docker','exec',name,'python','-c',reader],capture_output=True,text=True)
    print(run.stdout if run.returncode==0 else 'Container ended or log unavailable.',flush=True)
if not names:print('No running matching trial containers.')
