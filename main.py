from flask import Flask,jsonify,request
import diplom
from sympy import symbols,simplify
import json
from flask_cors import CORS, cross_origin
x = symbols('x')
app = Flask(__name__)
CORS(app)

@app.route("/")
def hello():
    return "Hello World! jslkjfsdlkj"

@app.route("/diplom", methods=['POST'])
def aproximation():
	data = json.loads(request.data)
	print(data)
	result = diplom.main(data['start'], data['end'], data['mu'],data['eps_mu'],simplify(data['func'] ))
	# print(result)
	return jsonify(result)
	# return jsonify(diplom.main(0.5,3.1,0.02,0.03,sin(x)))
	# return jsonify([[0.5, 1.8507812500000003, -0.137167050657332, 0.0122827611681035, 0.958761226399484, 0.0141201164446018], [1.8507812500000003, 2.9731262207031257, -0.0510442288070660, 0.0108070003159675, 0.147892293056749, 0.973927656225802], [2.9731262207031257, 3.1, -0.00581261089655455, 0.000113604222998573, -0.833699836503786, 2.79812191023127]])



if __name__ == "__main__":
    app.run(debug=True)
