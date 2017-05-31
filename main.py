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
    return ""

@app.route("/diplom", methods=['POST'])
def aproximation():
	data = json.loads(request.data)
	print(data)
	result = diplom.main(data['start'], data['end'], data['mu'],data['eps_mu'],simplify(data['func'] ),data['typ'])
	
	return jsonify(result)
	
if __name__ == "__main__":
    app.run(debug=True)
