// AABB only
// No angular momentum

"use strict";

class Vector {
	constructor(x /*: number */, y /*: number */) {
		this.x = x;
		this.y = y;
	}
}

function _vecDot(a /*: Vector */, b /*: Vector */) /*: number */ {
	return a.x * b.x + a.y * b.y;
}

function _vecPlus(a /*: Vector */, b /*: Vector */) /*: Vector */ {
	return new Vector(a.x + b.x, a.y + b.y);
}

function _vecMinus(a /*: Vector */, b /*: Vector */) /*: Vector */ {
	return new Vector(a.x - b.x, a.y - b.y);
}

function _vecMult(a /*: Vector */, s /*: number */) {
	return new Vector(a.x * s, a.y * s);
}

function _vecDiv(a /*: Vector */, s /*: number */) {
	return new Vector(a.x / s, a.y / s);
}

class AABB {
	constructor(min /*: Vector */, max /*: Vector */) {
		this._min = min;
		this._max = max;
	}
	
	x() {
		return this._min.x;
	}
	
	y() {
		return this._min.y;
	}
	
	width() {
		return this._max.x - this._min.x;
	}
	
	height() {
		return this._max.y - this._min.y;
	}
}

function _doCollide(a /*: AABB */, b /*: AABB */) /*: boolean */ {
	if (a._max.x < b._min.x || a._min.x > b._max.x) {
		return false;
	}
	if (a._max.y < b._min.y || a._min.y > b._max.y) {
		return false;
	}
	return true;
}

class Body {
	constructor(width /*: number */, height /*: number */) {
		this._aabb = new AABB(new Vector(0, 0), new Vector(width, height));
		this.velocity = new Vector(0, 0);
		this.position = new Vector(0, 0);
		this.inv_mass = 0;
	}
	
	aabb() {
		return new AABB(_vecPlus(this._aabb._min, this.position), _vecPlus(this._aabb._max, this.position));
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
		} else {
			this.inv_mass = 1 / m;
		}
		return this;
	}
	
	render(ctx /*: context */) {
		ctx.beginPath();
		const aabb /*: AABB */ = this.aabb();
		ctx.rect(aabb.x(), aabb.y(), aabb.width(), aabb.height());
		ctx.stroke();
	}
}

class Collision {
	constructor(a /*: Body */, b /*: Body */) {
		this.a = a;
		this.b = b;
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
	
	resolve() {
		const relativeVelocity /*: Vector */ = _vecMinus(this.b.velocity, this.a.velocity);
		const velocityAlongNormal /*: number */ = _vecDot(relativeVelocity, this._normal);
	
		if (velocityAlongNormal > 0) {
			return;
		}
	
		const e /*: number */ = 0.9;
	
		let j /*: number */ = -(1 + e) * velocityAlongNormal;
		j /= this.a.inv_mass + this.b.inv_mass;
	
		const impulse /*: Vector */ = _vecMult(this._normal, j);
		this.a.velocity = _vecMinus(this.a.velocity, _vecMult(impulse, this.a.inv_mass));
		this.b.velocity = _vecPlus(this.b.velocity, _vecMult(impulse, this.b.inv_mass));		
	}

	positionCorrection() {
		const percent /*: number */ = 0.2;
		const slop /*: number */ = 0.01;
		const correction /*: Vector */ = _vecMult(this._normal, Math.max(this._penetration - slop, 0) / (this.a.inv_mass + this.b.inv_mass) * percent);
		this.a.position = _vecMinus(this.a.position, _vecMult(correction, this.a.inv_mass));
		this.b.position = _vecPlus(this.b.position, _vecMult(correction, this.b.inv_mass));
	}
}

function _getSegmentOverlap(min1 /*: number */, max1 /*: number */, min2 /*: number */, max2 /*: number */) /*: number */ {
	return Math.max(0, Math.min(max1, max2) - Math.max(min1, min2))
}

function _getCollision(a /*: Body */, b /*: Body */) /*: Collision */ {
	if (a.inv_mass === 0 && b.inv_mass === 0) {
		return null;
	}
	
	let c /*: Collision */ = new Collision(a, b);
	
	const aAABB /*: AABB */ = a.aabb();
	const bAABB /*: AABB */ = b.aabb();
	
	const xOverlap /*: number */ = _getSegmentOverlap(aAABB._min.x, aAABB._max.x, bAABB._min.x, bAABB._max.x);
	if (xOverlap <= 0) {
		return null;
	}
	
	const yOverlap /*: number */ = _getSegmentOverlap(aAABB._min.y, aAABB._max.y, bAABB._min.y, bAABB._max.y);
	if (yOverlap <= 0) {
		return  null;
	}
		
	const n /*: Vector */ = _vecMinus(b.position, a.position);
			
	if (xOverlap < yOverlap) {
		if (n.x < 0) {
			c.setNormal(new Vector(-1, 0));
		} else {
			c.setNormal(new Vector(1, 0));
		}
		c.setPenetration(xOverlap);
	} else {
		if (n.y < 0) {
			c.setNormal(new Vector(0, -1));
		} else {
			c.setNormal(new Vector(0, 1));
		}
		c.setPenetration(yOverlap);
	}
	return c;
}

function step(bodies /*: Body[] */) {
	bodies
		.filter(b => b.inv_mass > 0)
		.forEach(b => {
			b.velocity = _vecPlus(b.velocity, new Vector(0, 1));
			b.position = _vecPlus(b.position, b.velocity);
		});
	
	for (let i = 0; i < bodies.length; ++i) {
		for (let j = i + 1; j < bodies.length; ++j) {
			const c /*: Collision */ = _getCollision(bodies[i], bodies[j]);
			if (c != null) {
				c.resolve();
				c.positionCorrection();
			}
		}
	}
}

function render(e /*: Canvas2d */, bodies /*: Body[] */) {
	const ctx = e.getContext('2d');
	ctx.clearRect(0, 0, e.width, e.height);
	bodies.forEach(b => b.render(ctx));
}