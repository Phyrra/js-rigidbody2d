// https://github.com/tutsplus/ImpulseEngine

"use strict";

class Vector {
	constructor(x /*: number */, y /*: number */) {
		this.x = x;
		this.y = y;
	}
	
	length() {
		return Math.sqrt(this.x * this.x + this.y * this.y);
	}
	
	normalize() {
		const length /*:number */ = this.length();
		this.x /= length;
		this.y /= length;
	}
}

function _vecDot(a /*: Vector */, b /*: Vector */) /*: number */ {
	return a.x * b.x + a.y * b.y;
}

function _vecCross(a /*: Vector */, b /*: Vector */) /*: number */ {
	return a.x * b.y - a.y * b.x;
}

//function _vecVectorCrossScalar(v /*: Vector */, s /*: number */) /*: Vector */ {
//	return new Vector(s * v.y, -s * v.x);
//}

function _vecScalarCrossVector(s /*: number */, v /*: Vector */) /*: Vector */ {
	return new Vector(-s * v.y, s * v.x);
}

function _vecPlus(a /*: Vector */, b /*: Vector */) /*: Vector */ {
	return new Vector(a.x + b.x, a.y + b.y);
}

function _vecMinus(a /*: Vector */, b /*: Vector */) /*: Vector */ {
	return new Vector(a.x - b.x, a.y - b.y);
}

function _vecMult(v /*: Vector */, s /*: number */) {
	return new Vector(v.x * s, v.y * s);
}

function _vecDiv(a /*: Vector */, s /*: number */) {
	return new Vector(v.x / s, v.y / s);
}

class Matrix {
	constructor() {
		this.m00 = 0; this.m01 = 0;
		this.m10 = 0; this.m11 = 0;
	}
	
	forRotation(r /*: number */) /*: Matrix */ {
		const c = Math.cos(r);
		const s = Math.sin(r);
		
		this.m00 = c; this.m01 = -s;
		this.m10 = s; this.m11 =  c;
		
		return this;
	}
	
	transpose() /*: Matrix */ {
		const m /*: Matrix */ = new Matrix();
		
		m.m00 = this.m00; m.m01 = this.m10;
		m.m10 = this.m01; m.m11 = this.m11;
		
		return m;
	}
}

function _matMult(m /*: Matrix */, v /*: Vector */) /*: Vector */ {
	return new Vector(m.m00 * v.x + m.m01 * v.y, m.m10 * v.x + m.m11 * v.y);
}

class Body {
	constructor(width /*: number */, height /*: number */) {
		this._width = width;
		this._height = height;
		this.velocity = new Vector(0, 0);
		this.position = new Vector(0, 0);
		this._center = new Vector(width / 2, height / 2);
		this.inv_mass = 0;
		this.inv_inertia = 0;
		this.rotation = 0;
		this.angularVelocity = 0;
		
		this._vertices /*: Vector[] */ = [
			new Vector(0, 0),
			new Vector(width, 0),
			new Vector(width, height),
			new Vector(0, height)
		];
		
		this._normals /*: Vector[] */ = [
			new Vector(0, -1),
			new Vector(1, 0),
			new Vector(0, 1),
			new Vector(-1, 0)
		];
		
		this._sprite = null;
	}
	
	vertices() /*: Vector[] */ {
		return this._vertices.map(v => _vecMinus(v, this._center));
	}
	
	normals() /*: Vector[] */ {
		return this._normals;
	}
	
	rotationMatrix() /*: Matrix */ {
		return new Matrix().forRotation(this.rotation);
	}
	
	setVelocity(v /*: Vector */) /*: Body */ {
		this.velocity = v;
		return this;
	}
	
	setPosition(p /*: Vector */) /*: Body */ {
		this.position = p;
		return this;
	}
	
	setMass(m /*: number */) /*: Body */ {
		if (m === 0) {
			this.inv_mass = 0;
			this.inv_inertia = 0;
		} else {
			this.inv_mass = 1 / m;
			this.inv_inertia = 1 / (m * this._width * this._height);
		}
		return this;
	}
	
	setCenter(c /*: Vector */) /*: Body */ {
		this._center = c;
		return this;
	}
	
	setRotation(r /*: number */) /*: Body */ {
		this.rotation = r;
		return this;
	}
	
	setAngularVelocity(w /*: number */) /*: Body */ {
		this.angularVelocity = w;
		return this;
	}
	
	setSprite(img /*: Image */) /*: Body */ {
		this._sprite = img;
		return this;
	}
	
	applyImpulse(impulse /*: Vector */, contactVector /*: Vector */) {
		if (this.inv_mass === 0) {
			return;
		}
		
		this.velocity = _vecPlus(this.velocity, _vecMult(impulse, this.inv_mass));
		this.angularVelocity = this.angularVelocity + this.inv_inertia * _vecCross(contactVector, impulse);
	}
	
	render(ctx /*: context */, offset /*: Vector */) {
		ctx.save();
		ctx.translate(offset.x, offset.y);
		ctx.translate(this.position.x, this.position.y);
		ctx.rotate(this.rotation);
		if (this._sprite == null) {
			ctx.beginPath();
			this.vertices().forEach(v => ctx.lineTo(v.x, v.y));
			ctx.closePath();
			ctx.stroke();
		} else {
			ctx.translate(-this._center.x, -this._center.y);
			ctx.drawImage(this._sprite, 0, 0, this._width, this._height);
		}
		ctx.restore();
	}
}

class Collision {
	constructor(a /*: Body */, b /*: Body */) {
		this.a = a;
		this.b = b;
		
		this._contactCount = 0;
		this._contacts = [
			new Vector(0, 0),
			new Vector(0, 0)
		];
		
		this._penetration = 0;
	}
	
	setNormal(n /*: Vector */) {
		this._normal = n;
	}
	
	normal() /*: Vector */ {
		return this._normal;
	}
	
	setPenetration(p /*: number */) {
		this._penetration = p;
	}
	
	penetration() /*: number */ {
		return this._penetration;
	}
	
	setContactCount(c /*: number */) {
		this._contactCount = c;
	}
	
	contactCount() /*: number */ {
		return this._contactCount;
	}
	
	setContact(i /*: number */, c /*: Vector */) {
		this._contacts[i] = c;
	}
	
	contact(i /*: number */) /*: Vector */ {
		return this._contacts[i];
	}
	
	resolve() {
		const restitution /*: number */ = 0.2;
		const staticFriction /*: number */ = 0.5;
		const dynamicFriction /*: number */ = 0.3;
	
		for (let i = 0; i < this._contactCount; ++i) {
			// Radius from center to contact
			const ra /*: Vector */ = _vecMinus(this._contacts[i], this.a.position);
			const rb /*: Vector */ = _vecMinus(this._contacts[i], this.b.position);
		
			//  Relative velocity
			let rv /*: Vector */ = _vecMinus(
				_vecPlus(this.b.velocity, _vecScalarCrossVector(this.b.angularVelocity, rb)),
				_vecMinus(this.a.velocity, _vecScalarCrossVector(this.a.angularVelocity, ra))
			);
		
			// Relative velocity along the normal
			const contactVel /*: number */ = _vecDot(rv, this._normal);
			if (contactVel > 0) {
				return;
			}
		
			const raCrossN = _vecCross(ra, this._normal);
			const rbCrossN = _vecCross(rb, this._normal);
			const invMassSum = this.a.inv_mass + this.b.inv_mass + raCrossN * raCrossN * this.a.inv_inertia + rbCrossN * rbCrossN * this.b.inv_inertia;
		
			// Impulse
			let j = -(1 + restitution) * contactVel;
			j /= invMassSum;
			j /= this._contactCount;
		
			// Apply impulse
			const impulse /*: Vector */ = _vecMult(this._normal, j);
			this.a.applyImpulse(_vecMult(impulse, -1), ra);
			this.b.applyImpulse(impulse, rb);
		
			// Friction impulse
			rv = _vecMinus(
				_vecPlus(this.b.velocity, _vecScalarCrossVector(this.b.angularVelocity, rb)),
				_vecMinus(this.a.velocity, _vecScalarCrossVector(this.a.angularVelocity, ra))
			);
		
			const t /*: Vector */ = _vecMinus(rv, _vecMult(this._normal, _vecDot(rv, this._normal)));
			t.normalize();
		
			// j tangent magnitude
			let jt = -_vecDot(rv, t);
			jt /= invMassSum;
			jt /= this._contactCount;
		
			if (Math.abs(jt) <= 1e-10) {
				return;
			}
		
			// Coulomb's law
			let tangentImpulse /*: Vector */;
			if (Math.abs(jt) < j * staticFriction) {
				tangentImpulse = _vecMult(t, jt);
			} else {
				tangentImpulse = _vecMult(t, -j * dynamicFriction);
			}
		
			// Apply tangent impulse
			this.a.applyImpulse(_vecMult(tangentImpulse, -1), ra);
			this.b.applyImpulse(tangentImpulse, rb);
		}
	}

	positionCorrection() {
		const slop /*: number */ = 0.05;
		const percent /*: number */ = 0.4;
		const correction /*: Vector */ = _vecMult(this._normal, Math.max(this._penetration - slop, 0) * percent / (this.a.inv_mass + this.b.inv_mass));
		this.a.position = _vecMinus(this.a.position, _vecMult(correction, this.a.inv_mass));
		this.b.position = _vecPlus(this.b.position, _vecMult(correction, this.b.inv_mass));
	}

	render(ctx /*: context */, offset /*: Vector */) {
		ctx.save();
		ctx.translate(offset.x, offset.y);
		if (this._contactCount === 0) {
			return;
		} else if (this._contactCount === 1) {
			ctx.beginPath();
			ctx.arc(this._contacts[0].x, this._contacts[0].y, 3, 0, 2 * Math.PI, false);
			ctx.fillStyle = '#ff0000';
			ctx.fill();
		} else {
			ctx.strokeStyle = '#ff0000';
			ctx.lineWidth = 3;
			ctx.beginPath();
			for (let i = 0; i < this._contactCount; ++i) {
				ctx.lineTo(this._contacts[i].x, this._contacts[i].y);
			}
			ctx.closePath();
			ctx.stroke();
			ctx.restore();
		}
		ctx.restore();
	}
}

function _getSupport(dir /*: Vector */, vertices /*: Vector[] */) /*: Vector */ {
	let bestProjection /*: number */ = Number.NEGATIVE_INFINITY;
	let bestVertex /*: Vector */;
	
	vertices.forEach(v => {
		const projection /*: number */ = _vecDot(v, dir);
		
		if (projection > bestProjection) {
			bestVertex = v;
			bestProjection = projection;
		}
	});
	
	return bestVertex;
}

function _findAxisLeastPenetration(a /*: Body */, b /*: Body */) /*: { number, number } */ {
	let bestDistance /*: number */ = Number.NEGATIVE_INFINITY;
	let bestIndex /*: number */ = -1;
	
	const aNormals /*: Vector[] */ = a.normals();
	const aVertices /*: Vector[] */ = a.vertices();
	const aRotation /*: Matrix */ = a.rotationMatrix();
	
	const bVertices /*: Vector[] */ = b.vertices();
	const bRotationT /*: Matrix */ = b.rotationMatrix().transpose();
	
	for (let i = 0; i < aNormals.length; ++i) {
		// Face normal of A
		let n /*: Vector */ = aNormals[i];
		n = _matMult(aRotation, n);
		
		// Transform face normal into B's model space
		n = _matMult(bRotationT, n);
		
		// Support point from B along -n
		const s /*: Vector */ = _getSupport(_vecMult(n, -1), bVertices);
		
		// Vertex on face from A, transform into B's model space
		let v /*: Vector */ = aVertices[i];
		v = _vecPlus(_matMult(aRotation, v), a.position);
		v = _vecMinus(v, b.position);
		v = _matMult(bRotationT, v);
		
		// Penetration distance (in B's model space)
		const d = _vecDot(n, _vecMinus(s, v));
		
		if (d > bestDistance) {
			bestDistance = d;
			bestIndex = i;
		}
	}
	
	return {
		index: bestIndex,
		distance: bestDistance
	};
}

function _findIncidentFace(reference /*: Body */, incident /*: Body */, referenceIndex /*: number */) /*: Vector[] */ {
	const incidentRotation /* Matrix */ = incident.rotationMatrix();
	
	let referenceNormal /*: Vector */ = reference.normals()[referenceIndex];
	
	// Normal in incident's frame of reference
	referenceNormal = _matMult(reference.rotationMatrix(), referenceNormal); // To world space
	referenceNormal = _matMult(incidentRotation.transpose(), referenceNormal); // To incident's model space
	
	// Find most anti-normal face on incident
	const incidentNormals /*: Vector[] */ = incident.normals();
	const incidentVertices /*: Vector[] */ = incident.vertices();

	let minDot /*: number */ = Number.POSITIVE_INFINITY;
	let incidentFace /*: number */ = -1;
	for (let i = 0; i < incidentNormals.length; ++i) {
		const dot =  _vecDot(referenceNormal, incidentNormals[i]);
		if (dot < minDot) {
			minDot = dot;
			incidentFace =  i;
		}
	}
	
	return [
		_vecPlus(_matMult(incidentRotation, incidentVertices[incidentFace]), incident.position),
		_vecPlus(_matMult(incidentRotation, incidentVertices[(incidentFace + 1) % incidentVertices.length]), incident.position)
	];
}

function _clip(n /*: Vector */, c /*: number */, face /*: Vector[] */) /*: { Vector[], number } */ {
	let sp /*: number */ = 0;
	const out /*: Vector[] */ = [
		face[0],
		face[1]
	];
	
	// d = ax + by - c
	const d1 /*: number */ = _vecDot(n, face[0]) - c;
	const d2 /*: number */ = _vecDot(n, face[1]) - c;
	
	// if negative (i.e. behind the plane), clip
	if (d1 <= 0) {
		out[sp++] = face[0];
	}
	if (d2 <= 0) {
		out[sp++] = face[1];
	}
	
	if (d1 * d2 < 0) {
		const alpha /*: number */ = d1 / (d1 - d2);
		out[sp] = _vecPlus(face[0], _vecMult(_vecMinus(face[1], face[0]), alpha));
		++sp;
	}
	
	if (sp === 3) {
		throw 'sp === 3';
	}
	
	return {
		face: [
			out[0],
			out[1]
		],
		count: sp
	};
}

function _biasGreaterThan(a /*: number */, b /*: number */) /*: boolean */ {
	const biasRelative /*:  number */ =  0.95;
	const biasAbsolute /*: number */ = 0.01;
	return a >= b * biasRelative + a * biasAbsolute;
}

function _getCollision(a /*: Body */, b /*: Body */) /*: Collision */ {
	if (a.inv_mass === 0 && b.inv_mass === 0) {
		return null;
	}
	
	const c /*: Collision */ = new Collision(a, b);
	c.setContactCount(0);
	
	const colA /*: { number, number } */ = _findAxisLeastPenetration(a, b);
	if (colA.distance >= 0) {
		return null;
	}
	
	const colB /*: { number, number } */ = _findAxisLeastPenetration(b, a);
	if (colB.distance >= 0) {
		return null;
	}
	
	let referenceIndex /*: number */;
	let flip /*: boolean */;
	
	let reference /*: Body */;
	let incident /*: Body */;
	if (_biasGreaterThan(colA.distance, colB.distance)) {
		reference = a;
		incident = b;
		referenceIndex = colA.index;
		flip = false;
	} else {
		reference = b;
		incident = a;
		referenceIndex = colB.index;
		flip = true;
	}
	
	let incidentFace = _findIncidentFace(reference, incident, referenceIndex);

	const referenceRotation /*: Matrix */ = reference.rotationMatrix();

	// Reference face vertices
	const referenceVertices /*: Vector[] */ = reference.vertices();	
	let v1 = referenceVertices[referenceIndex];
	let v2 = referenceVertices[(referenceIndex + 1) % referenceVertices.length];
	
	// To world
	v1 = _vecPlus(_matMult(referenceRotation, v1), reference.position);
	v2 = _vecPlus(_matMult(referenceRotation, v2), reference.position);
	
	const sidePlaneNormal /*: Vector */ = _vecMinus(v2, v1);
	sidePlaneNormal.normalize();
	
	const refFaceNormal /*: Vector */ = new Vector(sidePlaneNormal.y, -sidePlaneNormal.x);
	
	// ax + by = c
	const refC /*: number */ = _vecDot(refFaceNormal, v1);
	const negSide /*: number */ = -_vecDot(sidePlaneNormal, v1);
	const posSide /*: number */ = _vecDot(sidePlaneNormal, v2);
	
	const clip1 = _clip(_vecMult(sidePlaneNormal, -1), negSide, incidentFace);
	if (clip1.count < 2) {
		// Due to floating point error, possible to not have required points
		return null;
	}
	incidentFace = clip1.face;
	
	const clip2 = _clip(sidePlaneNormal, posSide, incidentFace);
	if (clip2.count <  2) {
		// Due to floating point error, possible to not have required points
		return null;
	}
	incidentFace = clip2.face;
	
	c.setNormal(flip ? _vecMult(refFaceNormal, -1) : refFaceNormal);
	
	// Keep points behind  reference face
	let cp /*: number */ = 0; // Clipped points behind reference face
	let separation /*: number */;
	
	separation = _vecDot(refFaceNormal, incidentFace[0]) - refC;
	if (separation <= 0) {
		c.setContact(cp, incidentFace[0]);
		c.setPenetration(-separation);
		++cp;
	}
	
	separation = _vecDot(refFaceNormal, incidentFace[1])  - refC;
	if (separation <= 0) {
		c.setContact(cp, incidentFace[1]);
		c.setPenetration(c.penetration() - separation);
		++cp;
	}
	
	c.setPenetration(c.penetration() / cp); // Average penetration
	c.setContactCount(cp);
	
	return c;
}

class World {
	constructor() {
		this._bodies = [];
		this._collisionIterations = 1;
		this._debug = false;
		this._hasGravity = false;
		this._gravity = new Vector(0, 0);
		this._offset = new Vector(0, 0);
	}
	
	setCollisionIterations(i /*: number */) /*: World */ {
		this._collisionIterations = i;
		return this;
	}
	
	setGravity(g /*: Vector */) /*: World */ {
		this._hasGravity = true;
		this._gravity = g;
		return this;
	}
	
	unsetGravity() /*: World */ {
		this._hasGravity = false;
		this._gravity = new Vector(0, 0);
		return this;
	}
	
	setBodies(bodies /*: Body[] */) /*: World */ {
		this._bodies = bodies.slice();
		return this;
	}
	
	addBody(b /*: Body */) /*: World */ {
		this._bodies.push(b);
		return this;
	}
	
	setDebug(d /*: boolean */) /*: World */ {
		this._debug = d;
		return this;
	}
	
	setOffset(o /*: Vector */) /*: World */ {
		this._offset = o;
		return this;
	}
	
	shiftBy(s /*: Vector */) /*: World */ {
		this._bodies.forEach(b => b.position = _vecMinus(b.position, s));
	}
	
	step() {
		const collisions /*: Collision[] */ = [];
		for (let i = 0; i < this._bodies.length; ++i) {
			for (let j = i + 1; j < this._bodies.length; ++j) {
				const c /*: Collision */ = _getCollision(this._bodies[i], this._bodies[j]);
				if (c != null) {
					collisions.push(c);
				}
			}
		}
		
		if (this._debug) {
			this._lastCollisions = collisions;
		}
	
		// Integrate forces
		if (this._hasGravity) {
			this._bodies
				.filter(b => b.inv_mass > 0)
				.forEach(b => b.velocity = _vecPlus(b.velocity, this._gravity));
		}
		// Resolve collisions
		for (let i = 0; i < this._collisionIterations; ++i) {
			collisions.forEach(c => c.resolve());
		}
		// Integrate velocity 
		this._bodies
			.filter(b => b.inv_mass > 0)
			.forEach(b => {
				b.position = _vecPlus(b.position, b.velocity);
				b.rotation = b.rotation + b.angularVelocity;
			});
		// Correct positions
		collisions.forEach(c => c.positionCorrection());
	}
	
	render(e /*: Canvas2d */) {
		const ctx = e.getContext('2d');
		ctx.clearRect(0, 0, e.width, e.height);
		this._bodies.forEach(b => b.render(ctx, this._offset));
		
		if (this._debug) {
			this._lastCollisions.forEach(c => c.render(ctx, this._offset));
		}
	}
}