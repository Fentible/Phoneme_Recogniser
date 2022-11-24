import os
import subprocess
import threading

from kivy.app import App
from kivy.lang.builder import Builder
from kivy.uix.screenmanager import ScreenManager, Screen
from kivy.uix.dropdown import DropDown

from enum import Enum

Builder.load_file('boxy.kv')

STOP = 1

# Main screen (options)
def stop_dtw():
    os.system('\n')

class Dataset(Enum):
    FEMALE: 'FEMALE'
    MALE: 'MALE'
    ALL: ''
    SPKR1: 'SPKR1'

class KNN(Enum):
    KNN: 'KNN'
    GROUP: 'GROUP'
    VOICED: 'VOICED'
    ZC: 'ZC'
    STE: 'STE'

class Bounds(Enum):
    BOUNDS: 'BOUNDS'
    NONE: ''
    
da = {
    'paa': '2', 'window': '128', 'banks': '40', 'paa_op': '0', 'dtw_window': '200',
    'interval_div': '2', 'nfft': '512', 'trunc': '24', 'mfccs': '99999', 'knn': '7',
    'group_k': '7', 'voice_k': '7', 'test_iter': '1', 
    }

bounds_args = {
    'zc_incr': '95', 'ste_incr': '85000', 'entr_incr': '400', 'neg_incr': '-300',
    'larg_incr': '1000',
    }


class MyGrid(Screen):

    def update_values(self):
        global da
        for i in self.ids:
            if i.endswith("in"):
                tmp = i[:-len("_in")]
                data = str(self.ids[str(i)].text)
                if data == "":
                    data = str(self.ids[str(i)].hint_text)
                # print(da[tmp])
                da[tmp] = data
                # print(da[tmp])

    def reset_values(self):
        global da
        for i in self.ids:
            if i.endswith("in"):
                tmp = i[:-len("_in")]
                self.ids[str(i)].text = ""
                self.ids[str(i)].hint_text = da[tmp]
                

    def change(self, *args):
        self.update_values()
        self.manager.current = "output"
        

# Output screen
class MyFrame(Screen):

    def runProcess(self):
        exe = ["dtw", 'paa', da['paa'], 'window', da['window'], 'banks', da['banks'], 'paa_op', da['paa_op'],
               'dtw_window', da['dtw_window'], 'interval_div', da['interval_div'], 'nfft', da['nfft'],
               'trunc', da['trunc'], 'mfccs', da['mfccs'], 'knn', da['knn'], 'group_k', da['group_k'],
               'voice_k', da['voice_k'], 'test_iter', da['test_iter']
               ]
        p = subprocess.Popen(['stdbuf', '-o0'] + exe, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        while True:
            if STOP == 1:
                p.terminate()
            retcode = p.poll()
            line = p.stdout.readline()
            self.update_out(line.decode("ascii"))
            # print(str(line))
            yield line
            if retcode is not None:
                break

    def update_out(self, text):
        self.ids.out.text += text

    def run_dtw(self):
        self.ids.out.text = ""
        for _ in self.runProcess():
            continue
        self.update_out("Finished :: Generating Matrix")
        subprocess.call(['python3', 'matrix.py'])
        self.ids.out.text = "Return to change values or press Start to run again..."
        self.ids.start_btn.text = "Start"

    def start_dtw_thread(self):
        
        global STOP
        if STOP == 1:
            STOP = 0
            self.ids.start_btn.text = "Stop"
            threading.Thread(target=self.run_dtw).start()
        else:
            STOP = 1
            self.ids.start_btn.text = "Start"

    def change(self, *args):
        self.manager.current = "main"


class MyApp(App):
    def build(self):
        sm = ScreenManager()
        main_screen = MyGrid(name='main')
        output_screen = MyFrame(name='output')
        sm.add_widget(main_screen)
        sm.add_widget(output_screen)
        return sm


if __name__ == "__main__":
    MyApp().run()
